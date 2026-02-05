#pragma once
#include "advpt/htgfile/VtkHtgFile.hpp"
#include "oktal/octree/CellGrid.hpp"
#include "oktal/octree/CellOctree.hpp"
#include "oktal/data/GridVector.hpp"
#include <filesystem>

namespace advpt::htgfile {
class SnapshotHtgFile;
}

namespace oktal::io::vtk {

void exportOctree(const CellOctree &octree,
                  const std::filesystem::path &filepath);

class CellGridExporter {
public:
  CellGridExporter(std::unique_ptr<advpt::htgfile::SnapshotHtgFile> htg_file,
                   const CellGrid &cells);

  template <typename T>
  CellGridExporter &writeGridVector(const std::string &name,
                                    std::span<const T> grid_vector);

  template <typename T, size_t Q>
  CellGridExporter &writeGridVector(const std::string &id,
                                    GridVectorView<const T, Q> vector);

private:
  std::unique_ptr<advpt::htgfile::SnapshotHtgFile> htg_file_;
  const CellGrid *cells_;
  std::size_t total_nodes_{0};
};

CellGridExporter exportCellGrid(const CellGrid &cells,
                                const std::filesystem::path &filepath);

template <typename T>
CellGridExporter &
CellGridExporter::writeGridVector(const std::string &name,
                                  std::span<const T> grid_vector) {
  std::vector<T> cell_data(total_nodes_, T{0});

  const auto morton_indices = cells_->mortonIndices();
  const auto &octree = cells_->octree();

  for (std::size_t enum_idx = 0; enum_idx < morton_indices.size(); ++enum_idx) {
    if (enum_idx < grid_vector.size()) {
      const MortonIndex morton_idx = morton_indices[enum_idx];
      const auto cell = octree.getCell(morton_idx);
      if (cell) {
        const std::size_t stream_idx = cell->streamIndex();
        cell_data[stream_idx] = grid_vector[enum_idx];
      }
    }
  }
  htg_file_->writeCellData<T>(name, cell_data);

  return *this;
}

template <typename T, size_t Q>
CellGridExporter &
CellGridExporter::writeGridVector(const std::string &id, GridVectorView<const T, Q> vector) {
  // Prepare cell data in AoS layout (layout_right)
  std::vector<T> cell_data(total_nodes_ * Q, T{0});

  const auto morton_indices = cells_->mortonIndices();
  const auto &octree = cells_->octree();

  for (std::size_t enum_idx = 0; enum_idx < morton_indices.size(); ++enum_idx) {
    if (enum_idx < vector.size()) {
      const MortonIndex morton_idx = morton_indices[enum_idx];
      const auto cell = octree.getCell(morton_idx);
      if (cell) {
        const std::size_t stream_idx = cell->streamIndex();

        // Copy Q components
        for (std::size_t q = 0; q < Q; ++q) {
          // SoA layout (grid_vector is layout_left)
          T value;
          if constexpr (Q == 1) {
            value = vector[enum_idx]; // 1D case
          } else {
            value = vector[enum_idx, q]; // 2D case
          }
          // AoS layout in cell_data
          cell_data[stream_idx * Q + q] = value;
        }
      }
    }
  }

  // Create mdspan view with AoS layout
  using MdCellDataExtents = std::experimental::dextents<std::size_t, 2>;
  using MdCellDataView = std::experimental::mdspan
                          <T,
                          MdCellDataExtents,
                          std::experimental::layout_right>;
                          
  const MdCellDataView md_cell_data_view(cell_data.data(), total_nodes_, Q);

  htg_file_->writeCellData<T>(id, md_cell_data_view);

  return *this;
}

} // namespace oktal::io::vtk
