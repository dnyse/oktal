/*
./poisson <refinementLevel> <max-iterations> <epsilon> <output-file>
*/

#include "oktal/io/VtkExport.hpp"
#include "oktal/octree/CellGrid.hpp"
#include "oktal/octree/CellOctree.hpp"
#include <cmath>
#include <filesystem>
#include <format>
#include <iostream>
#include <oktal/geometry/Vec.hpp>
#include <print>
#include <span>
#include <string>

namespace {
    
// Calculate the source term f at a given point
double calculate_f(const oktal::Vec3D& point) {
    // f(x,y,z) = 3*pi^2*cos(pi*x)*cos(pi*y)*cos(pi*z)
    return 3.0 * M_PI * M_PI * std::cos(M_PI * point[0]) * std::cos(M_PI * point[1]) * std::cos(M_PI * point[2]);
}

// Calculate the boundary condition u at a given point
double calculate_u(const oktal::Vec3D& point) {
    // u(x,y,z) = cos(pi*x)*cos(pi*y)*cos(pi*z)
    return std::cos(M_PI * point[0]) * std::cos(M_PI * point[1]) * std::cos(M_PI * point[2]);
}

// calculate the L2 norm of the residual
// ||r||_2 = (1 / N_cells) * sqrt(sum_i r_i^2)
double l2norm(const std::vector<double>& residual, std::size_t num_cells) {
    double sum = 0.0;
    for (const double val : residual) {
        sum += val * val;
    }
    return (1 / static_cast<double>(num_cells)) * std::sqrt(sum);
}

// Initialize source term f and boundary conditions u
void init_source_term_and_boundary_cond(
    const oktal::CellGrid &grid, std::vector<double> &f,
    const std::vector<oktal::Vec<std::ptrdiff_t, 3UL>> &neighbor_offsets,
    std::vector<bool> &is_boundary_vec, std::vector<double> &u,
    std::vector<double> &u_tmp) {
  for (const auto &cell : grid) {
    const oktal::Vec3D center = cell.center();
    const std::size_t enum_index = cell.enumerationIndex();

    f[enum_index] = calculate_f(center);

    // Check if cell is on boundary
    bool is_boundary = false;
    for (const auto &offset : neighbor_offsets) {
      if (!cell.neighbor(offset).isValid()) {
        is_boundary = true;
        break;
      }
    }

    // mark boundary status
    is_boundary_vec[enum_index] = is_boundary;

    if (is_boundary) {
      const double boundary_value = calculate_u(center);
      u[enum_index] = boundary_value;
      u_tmp[enum_index] = boundary_value;
    }
  }
}

// summation over neighbors
double sum_over_neighbors(
    const oktal::CellGrid::CellView &cell,
    const std::vector<oktal::Vec<std::ptrdiff_t, 3UL>> &neighbor_offsets,
    const std::vector<double> &u) {

    double sum_neighbors = 0.0;

    for (const auto &offset : neighbor_offsets) {
        const oktal::CellGrid::CellView neighbor_cell = cell.neighbor(offset);
        if (neighbor_cell.isValid()) {
            const std::size_t neighbor_enum_index =
                neighbor_cell.enumerationIndex();
            sum_neighbors += u[neighbor_enum_index];
        }
    }
  return sum_neighbors;
}

// Poisson solver using Jacobi iterative method
void poisson_solver(
    const double epsilon, const std::size_t maxIterations,
    const oktal::CellGrid &grid, std::vector<bool> &is_boundary_vec,
    const std::vector<oktal::Vec<std::ptrdiff_t, 3UL>> &neighbor_offsets,
    std::vector<double> &u, std::vector<double> &u_tmp, std::vector<double> &f,
    const double h, std::vector<double> &residual, const std::size_t num_cells,
    const std::filesystem::path &outputFile) {
  // Jacobi Solver
  std::size_t iteration = 0;
  double residual_norm = std::numeric_limits<double>::max();
  double sum_neighbors = 0.0;

  while (residual_norm > epsilon && iteration < maxIterations) {
    // for each cell in the grid
    for (const auto &cell : grid) {
      const std::size_t enum_index = cell.enumerationIndex();

      //  finite-difference stencil can only be applied at the interior cells
      if (is_boundary_vec[enum_index]) {
        continue; // skip boundary cells
      }

      sum_neighbors = sum_over_neighbors(cell, neighbor_offsets, u);

      // update u_tmp using Jacobi formula
      u_tmp[enum_index] = (f[enum_index] * h * h + sum_neighbors) / 6.0;
    }

    // swap u and u_tmp
    std::swap(u_tmp, u);

    // compute residual
    for (const auto &cell : grid) {
      const std::size_t enum_index = cell.enumerationIndex();

      if (is_boundary_vec[enum_index]) {
        continue;
      }

      sum_neighbors = sum_over_neighbors(cell, neighbor_offsets, u);

      residual[enum_index] =
          f[enum_index] - (6.0 * u[enum_index] - sum_neighbors) / (h * h);
    }

    iteration++;
    residual_norm = l2norm(residual, num_cells);
  }

  // Results output
  std::cout << "Jacobi solver completed in " << iteration
            << " iterations with residual norm " << residual_norm << "\n";

  // Write solution to output file
  oktal::io::vtk::exportCellGrid(grid, outputFile)
      .writeGridVector<double>("u", u)
      .writeGridVector<double>("f", f)
      .writeGridVector<double>("residual", residual);

  std::cout << "Solution successfully written to " << outputFile.string()
            << "\n";
}

} // namespace

int main(int argc, char* argv[]) {

    try {
        const std::span<char*> args(argv, static_cast<std::size_t>(argc));

        // input validation
        if (argc != 5) {
                std::println(std::clog, "Usage: {} <refinementLevel> <max-iterations> <epsilon> <output-file>", args[0]); 
                return 1;
        }

        // Parse command-line arguments
        const std::size_t refinementLevel = std::stoul(args[1]);
        const std::size_t maxIterations = std::stoul(args[2]);
        const double epsilon = std::stod(args[3]);
        const std::filesystem::path outputFile = args[4];

        // validate arguments
        if (refinementLevel < 0) {
            std::println(std::clog, "refinementLevel must be non-negative");
        }
        if (maxIterations <= 0) {
            std::println(std::clog, "max-iterations must be positive");
        }
        if (epsilon <= 0.0) {
            std::println(std::clog, "epsilon must be positive");
        }
        
        // Create octree grid
        const std::shared_ptr<const oktal::CellOctree> octree = oktal::CellOctree::createUniformGrid(refinementLevel);
        
        // define neighbor offsets for 6-connectivity
        const std::vector<oktal::Vec<std::ptrdiff_t, 3>> neighbor_offsets = {
            {-1, 0, 0}, // left, -X direction
            {1, 0, 0}, // right, +X direction
            {0, -1, 0}, // down, -Y direction
            {0, 1, 0}, // up, +Y direction
            {0, 0, -1}, // back, -Z direction
            {0, 0, 1} // front, +Z direction
        };

        // Create cell grid
        const oktal::CellGrid grid = oktal::CellGrid::create(octree)
                                    .levels({refinementLevel})
                                    .neighborhood(neighbor_offsets)
                                    .periodicityMapper(oktal::NoPeriodicity())
                                    .build();

        const std::size_t num_cells = grid.size(); // number of cells at the refinement level
        const std::size_t N = 1 << refinementLevel; // number of cells along one axis
        const double h = 1.0 / static_cast<double>(N); // cell spacing

        // set up grid vectors
        std::vector<double> f(num_cells, 0.0); // source term
        std::vector<double> u(num_cells, 0.0); // solution vector
        std::vector<double> u_tmp(num_cells, 0.0); // temporary solution
        std::vector<double> residual(num_cells, 0.0); // residual vector

        std::vector<bool> is_boundary_vec(num_cells, false); // is boundary or not

        // Initialize source term f and boundary conditions u
        init_source_term_and_boundary_cond(grid, f, neighbor_offsets,
                                           is_boundary_vec, u, u_tmp);

        poisson_solver(epsilon, maxIterations, grid, is_boundary_vec,
                       neighbor_offsets, u, u_tmp, f, h, residual, num_cells,
                       outputFile);

        return 0;

    } catch (const std::exception& e) {
        try {
            std::println(std::clog, "Error: {}", e.what());
        } catch (...) {
            return 1;
        }
        return 1;
    } catch (...) {
        try {
            std::println(std::clog, "An unknown error occurred.");
        } catch (...) {
            return 1;
        }
        return 1;
    }

    return 0;
}


