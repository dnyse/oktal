#include "oktal/octree/CellOctree.hpp"
#include <queue>
#include <stdexcept>

namespace oktal {
std::span<const CellOctree::Node> CellOctree::nodesStream() const {
  return std::span<const CellOctree::Node>{nodes_};
}
std::span<const CellOctree::Node>
CellOctree::nodesStream(std::size_t level) const {
  if (level >= level_start_idx_.size()) {
    return {};
  }
  const std::size_t start_index = level_start_idx_.at(level);
  const std::size_t count = nodes_per_level_.at(level);
  const std::span<const Node> complete_stream{nodes_};

  return complete_stream.subspan(start_index, count);
}

CellOctree::Node CellOctree::parseNodeChar(const char c) {
  bool refined = false;
  bool phantom = false;

  switch (c) {
  case '.': // Not Refined, Not Phantom
    break;
  case 'R': // Refined, Not Phantom
    refined = true;
    break;
  case 'P': // Not Refined, Phantom
    phantom = true;
    break;
  case 'X': // Refined, Phantom
    refined = true;
    phantom = true;
    break;
  default:
    throw std::invalid_argument("Invalid node symbol in descriptor.");
  }
  return {refined, phantom};
}

void CellOctree::validateSeperator(
    const std::size_t nodes_on_current_level, const std::size_t current_level,
    const std::size_t expected_nodes_current_level) {
  if (nodes_on_current_level == 0) {
    throw std::invalid_argument("Descriptor invalid: Empty level detected.");
  }
  // The root level (level 0) is always 1. All subsequent levels must be
  // a multiple of 8.
  if (current_level > 0 && nodes_on_current_level % 8 != 0) {
    throw std::invalid_argument(
        "Descriptor invalid: Level size not a multiple of 8 (after root).");
  }

  if (current_level > 0 &&
      nodes_on_current_level != expected_nodes_current_level) {
    throw std::invalid_argument(
        "Descriptor invalid: Incorrect number of children for previous "
        "level's refined nodes.");
  }
}

CellOctree CellOctree::fromDescriptor(std::string_view descriptor) {
  CellOctree octree;
  // Clear the default root node
  octree.nodes_.clear();
  octree.level_start_idx_.clear();
  octree.nodes_per_level_.clear();

  std::size_t current_level = 0;
  std::size_t nodes_on_current_level = 0;
  std::size_t expected_children_count = 0;
  std::size_t expected_nodes_current_level = 0; // Level 0 has no parent

  std::queue<std::size_t> refined_parents;

  octree.level_start_idx_.push_back(0);

  for (const char c : descriptor) {
    if (c == '|') {
      validateSeperator(nodes_on_current_level, current_level,
                        expected_nodes_current_level);

      // Finalize level data
      octree.nodes_per_level_.push_back(nodes_on_current_level);
      current_level++;
      nodes_on_current_level = 0;
      expected_nodes_current_level =
          expected_children_count; // Save for next level
      expected_children_count = 0;
      octree.level_start_idx_.push_back(octree.nodes_.size());

    } else {

      const CellOctree::Node node = parseNodeChar(c);

      if (current_level > 0 && (nodes_on_current_level % 8 == 0)) {

        const std::size_t children_start_idx = octree.nodes_.size();

        if (refined_parents.empty()) {
          throw std::invalid_argument("Descriptor invalid: Found children "
                                      "block without a refined parent");
        }

        // Get the parent index from the queue
        const std::size_t parent_idx = refined_parents.front();
        refined_parents.pop();

        // Check for over-refinement
        // parent must have been refined to have children
        if (!octree.nodes_.at(parent_idx).isRefined()) {
          throw std::invalid_argument("Descriptor invalid: Non-refined node "
                                      "expected, but children were added.");
        }

        // Set the childrenStartIndex for the parent node
        octree.nodes_.at(parent_idx).setChildrenStartIndex(children_start_idx);
      }

      octree.nodes_.push_back(node);

      if (node.isRefined()) {
        // If this new node is refined, it will be a parent in the next level
        refined_parents.push(octree.nodes_.size() - 1);
        expected_children_count += 8;
      }

      nodes_on_current_level++;
    }
  }
  if (!refined_parents.empty()) {
    throw std::invalid_argument("Descriptor invalid: Refined nodes found "
                                "without corresponding children block.");
  }
  // Validate last level
  if (current_level > 0 &&
      nodes_on_current_level != expected_nodes_current_level) {
    throw std::invalid_argument(
        "Descriptor invalid: Incorrect number of children for previous "
        "level's refined nodes.");
  }
  octree.nodes_per_level_.push_back(nodes_on_current_level);

  return octree;
}
std::optional<CellOctree::CellView>
CellOctree::getCell(MortonIndex m_idx) const {

  // Check if root cell
  if (m_idx.isRoot()) {
    return getRootCell();
  }

  std::size_t current_node_index = 0;
  const std::size_t target_level = m_idx.level();

  for (std::size_t level = 0; level < target_level; ++level) {

    const Node &parent_node = nodes_.at(current_node_index);

    if (!parent_node.isRefined()) {
      return std::nullopt;
    }

    const std::size_t branch_index = m_idx.getBranchIndex(level);

    current_node_index = parent_node.childIndex(branch_index);

    if (current_node_index >= nodes_.size()) {
      return std::nullopt;
    }
  }
  const Node &target_node = nodes_.at(current_node_index);
  if (target_node.isPhantom()) {
    return std::nullopt;
  }
  return CellView(target_node, m_idx, current_node_index, geometry());
}
bool CellOctree::cellExists(MortonIndex m_idx) const {
  return getCell(m_idx).has_value();
}
std::optional<CellOctree::CellView> CellOctree::getRootCell() const {
  if (nodes_.empty() || nodes_.at(0).isPhantom()) {
    return std::nullopt;
  }
  return CellOctree::CellView(nodes_.at(0), MortonIndex(), 0uz, geometry());
}

} // namespace oktal
