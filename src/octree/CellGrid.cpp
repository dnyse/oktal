#include "oktal/octree/CellGrid.hpp"
#include <algorithm>
#include <ranges>
#include <stdexcept>

namespace oktal {

namespace {

std::size_t findNeighborEnumIndex(
    const MortonIndex &cell_morton, const Vec<std::ptrdiff_t, 3> &offset,
    const std::function<std::optional<Vec<std::size_t, 3>>(Vec<std::ptrdiff_t, 3>,
                                                           std::size_t)>
        &periodicity_mapper,
    const CellOctree &octree, const CellGrid &grid) {
  const auto cell_coords = cell_morton.gridCoordinates();
  const std::size_t cell_level = cell_morton.level();

  Vec<std::ptrdiff_t, 3> neighbor_coords_signed;
  neighbor_coords_signed[0] =
      static_cast<std::ptrdiff_t>(cell_coords[0]) + offset[0];
  neighbor_coords_signed[1] =
      static_cast<std::ptrdiff_t>(cell_coords[1]) + offset[1];
  neighbor_coords_signed[2] =
      static_cast<std::ptrdiff_t>(cell_coords[2]) + offset[2];

  auto neighbor_coords_unsigned =
      periodicity_mapper(neighbor_coords_signed, cell_level);

  if (!neighbor_coords_unsigned.has_value()) {
    return CellGrid::NO_NEIGHBOR;
  }

  const MortonIndex neighbor_morton = MortonIndex::fromGridCoordinates(
      cell_level, neighbor_coords_unsigned.value());

  if (!octree.cellExists(neighbor_morton)) {
    return CellGrid::NO_NEIGHBOR;
  }

  auto neighbor_cell = octree.getCell(neighbor_morton);
  if (!neighbor_cell) {
    return CellGrid::NO_NEIGHBOR;
  }

  const std::size_t neighbor_enum_idx = grid.getEnumerationIndex(*neighbor_cell);
  if (neighbor_enum_idx == CellGrid::NOT_ENUMERATED) {
    return CellGrid::NO_NEIGHBOR;
  }

  return neighbor_enum_idx;
}

} // namespace

const CellOctree &CellGrid::octree() const { return *octree_; }

std::span<const MortonIndex> CellGrid::mortonIndices() const {
  return std::span<const MortonIndex>{morton_indices_};
}

std::size_t CellGrid::getEnumerationIndex(std::size_t streamIndex) const {
  auto it = stream_to_enum_.find(streamIndex);
  if (it != stream_to_enum_.end()) {
    return it->second;
  }
  return NOT_ENUMERATED;
}

std::size_t
CellGrid::getEnumerationIndex(const CellOctree::CellView &cell) const {
  return getEnumerationIndex(cell.streamIndex());
}

std::span<const std::size_t>
CellGrid::neighborIndices(Vec<std::ptrdiff_t, 3> offset) const {
  auto it = adjacency_lists_.find(offset);
  if (it == adjacency_lists_.end()) {
    throw std::out_of_range("Neighbor list not found for given offset");
  }
  return std::span<const std::size_t>{it->second};
}

CellGrid::CellView CellGrid::operator[](std::size_t index) const {
  return {index, this};
}

CellGrid::CellIterator CellGrid::begin() const { return {this, 0}; }

CellGrid::CellIterator CellGrid::end() const {
  return {this, morton_indices_.size()};
}

std::size_t CellGrid::size() const { return morton_indices_.size(); }

CellGrid::CellGridBuilder &
CellGrid::CellGridBuilder::levels(std::span<const std::size_t> lvls) {
  levels_.assign(lvls.begin(), lvls.end());
  return *this;
}

CellGrid::CellGridBuilder &
CellGrid::CellGridBuilder::levels(std::initializer_list<std::size_t> lvls) {
  levels_.assign(lvls.begin(), lvls.end());
  return *this;
}

CellGrid CellGrid::CellGridBuilder::build() const {
  CellGrid grid;
  grid.octree_ = octree_;

  std::vector<std::size_t> levels_to_enumerate;
  if (levels_.empty()) {
    for (std::size_t i = 0; i < octree_->numberOfLevels(); ++i) {
      levels_to_enumerate.push_back(i);
    }
  } else {
    levels_to_enumerate = levels_;
  }

  std::ranges::sort(levels_to_enumerate);

  std::size_t enum_index = 0;

  for (const std::size_t level : levels_to_enumerate) {
    for (auto cell : octree_->horizontalRange(level)) {
      grid.morton_indices_.push_back(cell.mortonIndex());
      grid.stream_to_enum_[cell.streamIndex()] = enum_index;
      enum_index++;
    }
  }

  if (!neighborhood_.empty()) {
    const std::size_t num_cells = grid.morton_indices_.size();

    for (const auto &offset : neighborhood_) {
      grid.adjacency_lists_[offset] =
          std::vector<std::size_t>(num_cells, NO_NEIGHBOR);
    }

    for (std::size_t i = 0; i < num_cells; ++i) {
      const MortonIndex &cell_morton = grid.morton_indices_[i];

      for (const auto &offset : neighborhood_) {
        const std::size_t neighbor_idx = findNeighborEnumIndex(
            cell_morton, offset, periodicity_mapper_, *octree_, grid);
        if (neighbor_idx != NO_NEIGHBOR) {
          grid.adjacency_lists_[offset][i] = neighbor_idx;
        }
      }
    }
  }

  return grid;
}

CellGrid::CellGridBuilder &CellGrid::CellGridBuilder::neighborhood(
    std::span<const Vec<std::ptrdiff_t, 3>> offsets) {
  neighborhood_.assign(offsets.begin(), offsets.end());
  return *this;
}

CellGrid::CellGridBuilder &CellGrid::CellGridBuilder::neighborhood(
    std::initializer_list<Vec<std::ptrdiff_t, 3>> offsets) {
  neighborhood_.assign(offsets.begin(), offsets.end());
  return *this;
}

Vec<double, 3> CellGrid::CellView::center() const {
  if (!isValid()) {
    throw std::logic_error("Cannot get center of invalid CellView");
  }
  const MortonIndex &morton = grid_->morton_indices_[enum_index_];
  return grid_->octree_->geometry().cellCenter(morton);
}

Box<double> CellGrid::CellView::boundingBox() const {
  if (!isValid()) {
    throw std::logic_error("Cannot get boundingBox of invalid CellView");
  }
  const MortonIndex &morton = grid_->morton_indices_[enum_index_];
  return grid_->octree_->geometry().cellBoundingBox(morton);
}

MortonIndex CellGrid::CellView::mortonIndex() const {
  if (!isValid()) {
    throw std::logic_error("Cannot get Morton index of invalid CellView");
  }
  return grid_->morton_indices_[enum_index_];
}

std::size_t CellGrid::CellView::level() const {
  if (!isValid()) {
    throw std::logic_error("Cannot get level of invalid CellView");
  }
  return grid_->morton_indices_[enum_index_].level();
}

CellGrid::CellView
CellGrid::CellView::neighbor(Vec<std::ptrdiff_t, 3> offset) const {
  if (!isValid()) {
    return {};
  }

  try {
    auto neighbors = grid_->neighborIndices(offset);
    const std::size_t neighbor_idx = neighbors[enum_index_];

    if (neighbor_idx == NO_NEIGHBOR) {
      return {};
    }

    return {neighbor_idx, grid_};
  } catch (const std::out_of_range &) {
    return {};
  }
}

} // namespace oktal
