#pragma once
#include "oktal/octree/MortonIndex.hpp"
#include "oktal/octree/OctreeGeometry.hpp"
#include <cstddef>
#include <span>
#include <string_view>
#include <vector>

namespace oktal {

class CellOctree {
public:
  // Nested classes are cool for hiding implementation details.
  class Node {
  private:
    // [1b refined, 1b phantom, 62b childrenindex] (childrenindex corresponds to
    // first child and the rest are stored in contiguos memory
    std::size_t data_;

    // Most significant bit (ULL = unsigned long long)
    // REFINED: All eight children are appended to the stream, ordered by their
    // Morton indices
    static constexpr std::size_t REFINED_MASK = 1ULL << 63;
    // Second most significant bit
    // PHANTOM: Cells that are not required by the simulation
    static constexpr std::size_t PHANTOM_MASK = 1ULL << 62;
    // Mask for the index value
    static constexpr std::size_t INDEX_MASK = ~(REFINED_MASK | PHANTOM_MASK);

  public:
    // Default constructor: non-refined, non-phantom node.
    Node() : data_(0) {}

    // Constructor with all attributes. childrenIdx  0 as it's undefined for
    // non-refined.
    Node(bool refined, bool phantom, std::size_t childrenIdx = 0)
        : data_(childrenIdx & INDEX_MASK) {
      if (refined) {
        data_ |= REFINED_MASK;
      }
      if (phantom) {
        data_ |= PHANTOM_MASK;
      }
    }
    [[nodiscard]] bool isRefined() const { return (data_ & REFINED_MASK) != 0; }
    [[nodiscard]] bool isPhantom() const { return (data_ & PHANTOM_MASK) != 0; }
    [[nodiscard]] std::size_t childrenStartIndex() const {
      return data_ & INDEX_MASK;
    }

    void setRefined(bool refined) {
      if (refined) {
        data_ |= REFINED_MASK;
      } else {
        data_ &= ~REFINED_MASK;
      }
    }
    void setPhantom(bool phantom) {
      if (phantom) {
        data_ |= PHANTOM_MASK;
      } else {
        data_ &= ~PHANTOM_MASK;
      }
    }
    void setChildrenStartIndex(std::size_t index) {
      data_ = (data_ & ~INDEX_MASK) | (index & INDEX_MASK);
    }

    // Calculates the stream index for a child node given its octal branch index
    // (0-7).
    [[nodiscard]] std::size_t childIndex(size_t branch) const {
      // Children are always grouped as 8 siblings, starting at
      // childrenStartIndex()
      return childrenStartIndex() + branch;
    }
  };

  class CellView {
  private:
    Node node_;
    MortonIndex morton_index_;
    std::size_t stream_index_;
    OctreeGeometry geometry_;

  public:
    CellView() = default;
    CellView(Node node, MortonIndex m_idx, std::size_t stream_idx,
             const OctreeGeometry &geom)
        : node_(node), morton_index_(m_idx), stream_index_(stream_idx),
          geometry_(geom) {}
    [[nodiscard]] MortonIndex mortonIndex() const { return morton_index_; }
    [[nodiscard]] bool isRoot() const { return morton_index_.isRoot(); }
    [[nodiscard]] bool isRefined() const { return node_.isRefined(); }
    [[nodiscard]] size_t level() const { return morton_index_.level(); }
    [[nodiscard]] size_t streamIndex() const { return stream_index_; }
    [[nodiscard]] Vec<double, 3> center() const {
      return geometry_.cellCenter(morton_index_);
    }
    [[nodiscard]] Box<double> boundingBox() const {
      return geometry_.cellBoundingBox(morton_index_);
    }
  };

  class OctreeCursor {
  public:
    // default constructor
    OctreeCursor() : octree_(nullptr) {}

    // constructor with octree reference
    OctreeCursor(const CellOctree& octree) : octree_(&octree), path_{0} {}
    
    // constructor with octree reference and path
    OctreeCursor(const CellOctree& octree, std::span<const std::size_t> path) : octree_(&octree) {
      if (!path.empty()) {
        path_.assign(path.begin(), path.end());
      }
    }
     
    [[nodiscard]] const CellOctree* octree() const {
      return octree_;
    }

    [[nodiscard]] std::span<const std::size_t> path() const {
      return {path_.data(), path_.size()};
    }

    // Observers

    [[nodiscard]] bool empty() const {
      return octree_ == nullptr;
    }

    [[nodiscard]] bool end() const {
      return octree_ != nullptr && path_.empty();
    }

    [[nodiscard]] std::size_t currentLevel() const {
      return path_.empty() ? 0 : path_.size() - 1;
    }

    [[nodiscard]] std::size_t currentStreamIndex() const {
      return path_.back();
    }

    [[nodiscard]] std::optional<CellOctree::CellView> currentCell() const {
      if (end() || empty()) {
        return std::nullopt;
      }
      return octree_->getCell(mortonIndex());
    }

    [[nodiscard]] bool firstSibling() const {
      if (end() || empty()) {
        return false;
      }
      return mortonIndex().isFirstSibling();
    }

    [[nodiscard]] bool lastSibling() const {
      if (end() || empty()) {
        return false;
      }
      return mortonIndex().isLastSibling();
    }

    [[nodiscard]] MortonIndex mortonIndex() const {
      if (end() || empty()) {
        return {};
      }

      std::vector<morton_bits_t> branches;
      branches.reserve(path_.size());

      for (std::size_t i = 1; i < path_.size(); ++i) {
        const Node& parent = octree_->nodesStream()[path_[i - 1]]; 
        std::size_t branch = path_[i] - parent.childrenStartIndex();

        // check branch is valid
        if (branch > 7) {
          return {};
        }

        branches.push_back(static_cast<morton_bits_t>(branch));
      }

      return MortonIndex::fromPath(branches);
    }

    [[nodiscard]] bool operator==(const OctreeCursor& other) const {
      return octree_ == other.octree_ && path_ == other.path_;
    }

    [[nodiscard]] bool operator!=(const OctreeCursor& other) const {
      return !(*this == other);
    }

    // Traversal

    void ascend() {
      if (!path_.empty()) {
        path_.pop_back();
      }
    }

    void descend() {
      if (end() || empty()) {
        return;
      }

      const Node& currentNode = octree_->nodesStream()[currentStreamIndex()];
      if (currentNode.isRefined()) {
        path_.push_back(currentNode.childrenStartIndex());
      }      
    }

    void descend(std::size_t childIdx) {
      if (childIdx > 7) {
        throw std::out_of_range("Child index must be between 0 and 7.");
      }

      if (end() || empty()) {
        return;
      }

      const Node& currentNode = octree_->nodesStream()[currentStreamIndex()];
      if (currentNode.isRefined()) {
        path_.push_back(currentNode.childIndex(childIdx));
      }      
    }

    void previousSibling() {
      if (end() || empty() || mortonIndex().isFirstSibling()) {
        return;
      }
      
      const std::size_t currentIdx = currentStreamIndex();
      const Node& parentNode = octree_->nodesStream()[path_[path_.size() - 2]];
      const std::size_t firstSiblingIdx = parentNode.childIndex(0);
    
      if (currentIdx > firstSiblingIdx) {
        path_.back() -= 1;
      }
    }

    void nextSibling() {
      if (end() || empty() || mortonIndex().isLastSibling()) {
        return;
      }
      
      const std::size_t currentIdx = currentStreamIndex();
      const Node& parentNode = octree_->nodesStream()[path_[path_.size() - 2]];
      const std::size_t lastSiblingIdx = parentNode.childIndex(7);

      if (currentIdx < lastSiblingIdx) {
        path_.back() += 1;
      }
    }

    void toSibling(std::size_t siblingIdx) {
      if (end() || empty()) {
        return;
      }
      if (siblingIdx > 7) {
        throw std::out_of_range("Sibling index must be between 0 and 7.");
      }
     
      if (path_.size() == 1) {    // Root node
        if (siblingIdx != 0) {
          throw std::out_of_range("Root node has no siblings.");
        }

        return;
      }

      const Node& parentNode = octree_->nodesStream()[path_[path_.size() - 2]];
      path_.back() = parentNode.childIndex(siblingIdx);
      
    }

    void toEnd() {
      path_.clear();
    }

  private:
    const CellOctree* octree_;
    std::vector<std::size_t> path_;
  };

private:
  std::vector<Node> nodes_ = {{}};
  std::vector<std::size_t> level_start_idx_ = {0};
  std::vector<std::size_t> nodes_per_level_ = {1};
  OctreeGeometry geometry_;

public:
  CellOctree() = default;

  CellOctree(OctreeGeometry geometry) : geometry_(geometry) {}

  [[nodiscard]] std::size_t numberOfNodes() const { return nodes_.size(); }
  [[nodiscard]] std::size_t numberOfLevels() const {
    return nodes_per_level_.size();
  }
  [[nodiscard]] std::size_t numberOfNodes(std::size_t level) const {

    if (level >= nodes_per_level_.size()) {
      return 0; // TODO: Throw an error?
    }
    return nodes_per_level_.at(level);
  }
  [[nodiscard]] std::span<const Node> nodesStream() const;
  [[nodiscard]] std::span<const Node> nodesStream(std::size_t level) const;

  static CellOctree fromDescriptor(std::string_view descriptor);

  static Node parseNodeChar(char c);
  static void validateSeperator(std::size_t nodes_on_current_level,
                                std::size_t current_level,
                                std::size_t expected_nodes_current_level);
  [[nodiscard]] const OctreeGeometry &geometry() const { return geometry_; }
  [[nodiscard]] std::optional<CellView> getCell(MortonIndex m_idx) const;
  [[nodiscard]] bool cellExists(MortonIndex m_idx) const;
  [[nodiscard]] std::optional<CellView> getRootCell() const;

  [[nodiscard]] auto preOrderDepthFirstRange() const;
  [[nodiscard]] auto horizontalRange(std::size_t level) const;

};

// alias
using OctreeCursor = CellOctree::OctreeCursor;

// concept OctreeIteratorPolicy
template<typename T>
concept OctreeIteratorPolicy = 
  std::semiregular<T> &&
  requires(T const policy, OctreeCursor& cursor)
  {
    { policy.advance(cursor) } -> std::same_as<void>;
  };  

// iterator class template
template <OctreeIteratorPolicy TPolicy>
class OctreeIterator {
public:
  // constructor
  OctreeIterator() = default;

  OctreeIterator(OctreeCursor cursor, TPolicy policy)
      : cursor_(std::move(cursor)), policy_(std::move(policy)) {}

  // dereference operator
  [[nodiscard]] CellOctree::CellView operator*() const {
    auto cellOpt = cursor_.currentCell();
    if (!cellOpt.has_value()) {
      throw std::logic_error("Dereferencing an invalid iterator.");
    }
    return *cellOpt;
  }

  // pre-increment operator
  OctreeIterator& operator++() {
    policy_.advance(cursor_);
    return *this;
  }

  // post-increment operator
  OctreeIterator operator++(int) {
    OctreeIterator temp = *this;
    ++(*this);
    return temp;
  }

  // equality operator
  [[nodiscard]] bool operator==(const OctreeIterator& other) const {
    return cursor_ == other.cursor_;
  }

  // inequality operator
  [[nodiscard]] bool operator!=(const OctreeIterator& other) const {
    return !(*this == other);
  }

  // alias for iterator
  using iterator_category = std::forward_iterator_tag;
  using value_type = CellOctree::CellView;
  using difference_type = std::ptrdiff_t;
  using pointer = const CellOctree::CellView*;
  using reference = const CellOctree::CellView&;  

private:
    OctreeCursor cursor_;
    TPolicy policy_;
};

// range class template
template <OctreeIteratorPolicy TPolicy>
class OctreeCellsRange {
public:
  // constructor
  OctreeCellsRange() = default;

  OctreeCellsRange(OctreeCursor start, OctreeCursor end, TPolicy policy)
      : start_(std::move(start)), end_(std::move(end)), policy_(std::move(policy)) {}

  // begin iterator
  [[nodiscard]] OctreeIterator<TPolicy> begin() const {
    return OctreeIterator<TPolicy>(start_, policy_);
  }

  // end iterator
  [[nodiscard]] OctreeIterator<TPolicy> end() const {
    return OctreeIterator<TPolicy>(end_, policy_);
  }

private:
  OctreeCursor start_;
  OctreeCursor end_;
  TPolicy policy_;
};

// DFS Policy
class PreOrderDepthFirstPolicy {
public:

  void advance(OctreeCursor& cursor) const {

    if (cursor.end() || cursor.empty()) {
      return;
    }

    advanceToNextValid(cursor);

    while (!cursor.end() && !cursor.currentCell().has_value()) {
      advanceToNextValid(cursor);
    }

  }

private:
  // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
  void advanceToNextValid(OctreeCursor& cursor) const {

    if (cursor.end() || cursor.empty()) {
      return;
    }

    const auto& currentNode = cursor.octree()->nodesStream()[cursor.currentStreamIndex()];

    if (currentNode.isRefined()) {
      cursor.descend();
      return;
    }

    while (true) {
      if (!cursor.lastSibling()) {
        cursor.nextSibling();
        break;
      }
      
      cursor.ascend(); // Go up to parent

      if (cursor.end() || cursor.empty()) {
        return;
      }
    }
  }

};

inline auto CellOctree::preOrderDepthFirstRange() const {

  OctreeCursor start{*this}; // root cursor
  const PreOrderDepthFirstPolicy policy;
  
  while (!start.end() && !start.currentCell().has_value()) {
    policy.advance(start);
  }

  return OctreeCellsRange<PreOrderDepthFirstPolicy>(
    start,
    OctreeCursor(*this, {}), // end cursor
    policy
  );

}

// Horizontal Policy
class HorizontalPolicy {
public:
  HorizontalPolicy() = default;
  HorizontalPolicy(std::size_t level) : target_level_(level) {}

  void advance(OctreeCursor& cursor) const {
    if (cursor.end() || cursor.empty()) {
      return;
    }

    do {
      advanceToNextValid(cursor);
    } while (!cursor.end() && !cursor.currentCell().has_value());
  }

private:
  std::size_t target_level_;

  void advanceToNextValid(OctreeCursor& cursor) const {
    while (true) {
      if (!cursor.lastSibling()) {
        cursor.nextSibling();

        // dive down to target level
        if (diveToLevel(cursor)) {
          return;
        }
      } 
      else {
        cursor.ascend(); // Go up to parent
        if (cursor.path().empty()) {
          cursor.toEnd();
          return;
        }
      }
    }
  }

  [[nodiscard]] bool diveToLevel(OctreeCursor& cursor) const {
    while (cursor.currentLevel() <  target_level_) {
      const auto& node = cursor.octree()->nodesStream()[cursor.currentStreamIndex()];
      if (node.isRefined()) {
        cursor.descend(0); 
      } else {
        return false;
      }
    }
    return true; 
  }
};

inline auto CellOctree::horizontalRange(std::size_t level) const {
  OctreeCursor start(*this); 
  const HorizontalPolicy policy(level);

  if (level > 0) {
    bool perfect_path = true;
    OctreeCursor temp = start; 
    
    for (std::size_t i = 0; i < level; ++i) {
       const auto& node = nodesStream()[temp.currentStreamIndex()];
       if (node.isRefined()) {
         temp.descend(0);
       } else {
         perfect_path = false; 
         break; 
       }
    }

    start = temp;

    // if target level not reached, advance
    if (!perfect_path) {
       policy.advance(start);
    }
  }

  // advance to first valid cell
  while (!start.end() && !start.currentCell().has_value()) {
      policy.advance(start);
  }

  return OctreeCellsRange<HorizontalPolicy>(
      start,
      OctreeCursor(*this, {}), // end cursor
      policy
  );
}

} // namespace oktal
