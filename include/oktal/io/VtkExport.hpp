#pragma once
#include "oktal/octree/CellOctree.hpp"
#include <filesystem>

namespace oktal::io::vtk {

void exportOctree(const CellOctree &octree,
                  const std::filesystem::path &filepath);

} // namespace oktal::io::vtk
