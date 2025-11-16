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
  [[nodiscard]] const OctreeGeometry& geometry() const { return geometry_; }
  [[nodiscard]] std::optional<CellView> getCell(MortonIndex m_idx) const;
  [[nodiscard]] bool cellExists(MortonIndex m_idx) const;
  [[nodiscard]] std::optional<CellView> getRootCell() const;
};

} // namespace oktal
