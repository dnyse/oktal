#include "oktal/io/VtkExport.hpp"
#include "advpt/htgfile/VtkHtgFile.hpp"
#include "oktal/geometry/Box.hpp"

namespace oktal::io::vtk {
using advpt::htgfile::HyperTree;

// WARN: Because the bit streams are packed (for both mask and description) in
// uint8_t including the root node, the byte boundary falls between the 7th and
// 8th child. So the last child can always be found in the next uint8_t
const auto pack_bits =
    [](const std::vector<uint8_t> &bits) -> std::vector<uint8_t> {
  std::vector<uint8_t> bytes;
  if (bits.empty()) {
    return bytes;
  }
  const size_t num_bytes = (bits.size() + 7) / 8;
  bytes.resize(num_bytes, 0);

  for (size_t i = 0; i < bits.size(); ++i) {
    if (bits[i] == 1) {
      const size_t bit_index = i;
      const size_t byte_index = bit_index / 8;
      const size_t bit_offset = 7 - (bit_index % 8);
      bytes[byte_index] |= (1u << bit_offset);
    }
  }
  return bytes;
};

const auto createHyperTree = [](const CellOctree &octree,
                                std::vector<int> &level_data) -> HyperTree {
  HyperTree hyperTree;
  const auto &geometry = octree.geometry();
  const auto root_box = geometry.cellBoundingBox(MortonIndex());
  const Box<double>::vector_type min_corner = root_box.minCorner();
  const Box<double>::vector_type max_corner = root_box.maxCorner();
  hyperTree.xCoords = {min_corner[0], max_corner[0]};
  hyperTree.yCoords = {min_corner[1], max_corner[1]};
  hyperTree.zCoords = {min_corner[2], max_corner[2]};
  std::vector<uint8_t> descriptor_bits;
  std::vector<uint8_t> mask_bits;
  std::vector<MortonIndex> pseudo_phantom;

  const std::size_t levels = octree.numberOfLevels();
  for (std::size_t level = 0; level < levels; ++level) {
    const auto level_nodes = octree.nodesStream(level);
    hyperTree.nodesPerDepth.emplace_back(level_nodes.size());
    for (const auto node : level_nodes) {
      if (level < levels - 1) {
        descriptor_bits.push_back(node.isRefined() ? 1 : 0);
      }
      mask_bits.push_back(node.isPhantom() && !node.isRefined() ? 1 : 0);
      level_data.push_back(static_cast<int>(level));
    }
  }
  hyperTree.descriptor = pack_bits(descriptor_bits);
  hyperTree.mask = pack_bits(mask_bits);

  return hyperTree;
};
void exportOctree(const CellOctree &octree,
                  const std::filesystem::path &filepath) {
  std::vector<int> level_data;
  const HyperTree hyperTree = createHyperTree(octree, level_data);

  auto htg_file = advpt::htgfile::SnapshotHtgFile::create(filepath, hyperTree);
  htg_file.writeCellData<int>("level", level_data);
}

} // namespace oktal::io::vtk
