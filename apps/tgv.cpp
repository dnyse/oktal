/*
Taylor-Green Vortex simulation using Lattice Boltzmann Method (LBM).
./tgv <refinement-level> <output-directory>
*/

#include <iostream>
#include <span>
#include <filesystem>
#include <print>
#include <string>
#include <chrono>
#include <fstream>
#include "oktal/io/VtkExport.hpp"
#include "oktal/octree/CellGrid.hpp"
#include "oktal/octree/CellOctree.hpp"
#include <oktal/geometry/Vec.hpp>
#include "oktal/lbm/D3Q19.hpp"
#include "oktal/lbm/TaylorGreen.hpp"
#include "oktal/lbm/LbmKernels.hpp"

int main(int argc, char* argv[]) {
    try {
        const std::span<char*> args(argv, static_cast<std::size_t>(argc));

        // input validation
         if (argc != 3) {
            std::println(std::clog, "Usage: {} <refinement-level> <output-directory>", args[0]); 
            return 1;
         }

        // Parse command-line arguments
        const std::size_t refinementLevel = std::stoul(args[1]);
        const std::filesystem::path outputFile = args[2];

        // validate arguments
        // refinement level must be at least five
        if (refinementLevel < 5) {
            std::println(std::clog, "refinement-level must be at least 5");
            return 1;
        }

        // create output directory
        std::filesystem::create_directories(outputFile);

        // TGV simulation
        
        // Create octree
        const std::shared_ptr<const oktal::CellOctree> octree = oktal::CellOctree::createUniformGrid(refinementLevel);
        
        // Create cell grid
        const oktal::CellGrid grid = oktal::CellGrid::create(octree)
                                    .levels({refinementLevel})
                                    .neighborhood(oktal::lbm::D3Q19::CS_NO_CENTER)
                                    .periodicityMapper(oktal::Torus({true,true,true})) // fully periodic in all directions
                                    .build();

        // create LBM fields
        oktal::lbm::D3Q19Lattice f(grid); // current pdfs
        oktal::lbm::D3Q19Lattice fTmp(grid); // temporary pdfs
        oktal::GridVector<double,1> rho(grid); // density field
        oktal::GridVector<double,3> u(grid); // velocity field

        std::println(std::clog, "Grid setup complete:");
        std::println(std::clog, "  Refinement level: {}", refinementLevel);
        std::println(std::clog, "  Number of cells: {}", grid.size());
        std::println(std::clog, "  Output directory: {}", outputFile.string());

        // initial flow state

        // create an instance of TaylorGreen
        constexpr double latticeViscosity = 1.0 / 6.0;
        const oktal::lbm::TaylorGreen tgv(refinementLevel, latticeViscosity);

        // initialize density and velocity fields at t=0
        for (const auto& cell : grid) {
            const std::size_t cellIdx{cell};
            const auto center_position = cell.center();
            rho[cellIdx, 0] = tgv.rho(center_position * 2 * M_PI, 0);
            const auto velocity = tgv.u(center_position * 2 * M_PI, 0);
            u[cellIdx, 0] = velocity[0];
            u[cellIdx, 1] = velocity[1];
            u[cellIdx, 2] = velocity[2];
        }

        // initialize pdfs
        oktal::lbm::InitializePdfs initPdfs;
        for (const auto& cell : grid) {
            initPdfs(cell, f.view(), rho.const_view(), u.const_view());
        }

        std::println(std::clog, "flow state initialized at t=0");

        // integration loop
        const double omega = tgv.omega();
        
        oktal::lbm::ComputeMacroscopicQuantities computeMacroscopic;
        const oktal::lbm::Collide collide(omega);
        oktal::lbm::Stream stream;

        // run simulation
        std::println(std::clog, "Starting LBM loop:");
        const std::size_t numTimesteps = tgv.numberOfTimesteps();
        std::println(std::clog, "  Number of timesteps: {}", numTimesteps);
        constexpr std::size_t interval = 50; // export every 50 timesteps

        // start time
        const auto startTime = std::chrono::steady_clock::now();
        auto lastTime = startTime;
        
        for (std::size_t t = 0; t < numTimesteps; ++t) {
            // compute macroscopic quantities
            for (const auto& cell : grid) {
                computeMacroscopic(cell, f.const_view(), rho.view(), u.view());
            }

            // export intermediate results
            if (t % interval == 0) {
                const std::filesystem::path timeSeriesPath = outputFile / 
                    std::format("step{}.vtkhdf", t / interval);

                auto exporter = oktal::io::vtk::exportCellGrid(grid, timeSeriesPath);
                exporter.writeGridVector<double, 1>("density", rho.const_view());
                exporter.writeGridVector<double, 3>("velocity", u.const_view());

                std::println(std::clog, "  Exported results at timestep {} to {}", t, timeSeriesPath.string());
            }

            // collide
            for (const auto& cell : grid) {
                collide(cell, f.view(), rho.const_view(), u.const_view());            
            }

            // stream
            for (const auto& cell : grid) {
                stream(cell, fTmp.view(), f.const_view());
            }

            // swap
            std::swap(f, fTmp);

            // output every 10% of the simulation
            if ((t + 1) % interval == 0 || t == 0) {
                const auto currentTime = std::chrono::steady_clock::now();
                const auto elapsedTime = std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime).count();
                const auto intervalTime = std::chrono::duration_cast<std::chrono::seconds>(currentTime - lastTime).count();
                const double avgTimePerStep = static_cast<double>(elapsedTime) / static_cast<double>(t + 1);
                const std::size_t stepsRemaining = numTimesteps - (t + 1);
                const auto estimatedRemainingTime = static_cast<std::size_t>(avgTimePerStep * static_cast<double>(stepsRemaining));
                std::println(std::clog, "  Progress: {}/{} steps ({:.1f}%) | Elapsed: {:.2f}s | Interval: {:.2f}s | ETA: {:.2f}s",
                            t + 1,
                            numTimesteps, 
                            100.0 * static_cast<double>(t + 1) / static_cast<double>(numTimesteps),
                            static_cast<double>(elapsedTime) / 1000.0,
                            static_cast<double>(intervalTime) / 1000.0,
                            static_cast<double>(estimatedRemainingTime) / 1000.0);
                lastTime = currentTime;
            }
        }

        // end time
        const auto endTime = std::chrono::steady_clock::now();
        const auto totalElapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
        std::println(std::clog, "Simulation complete.");
        std::println(std::clog, "Total simulation time: {:.3f}s", static_cast<double>(totalElapsedTime.count()) / 1000.0);
    
        // export results to VTK
        // compute final macroscopic quantities
        for (const auto& cell : grid) {
            computeMacroscopic(cell, f.const_view(), rho.view(), u.view());
        }

        // compute time from timesteps
        const double finalTime = static_cast<double>(numTimesteps) * tgv.dt();

        // compute velocity errors
        oktal::GridVector<double,3> velocityError(grid);
        double error_L2_ux = 0.0;
        double error_L2_uy = 0.0;
        double norm_L2_ux = 0.0;
        double norm_L2_uy = 0.0;

        for (const auto& cell : grid) {
            const std::size_t cellIdx{cell};
            const auto center_position = cell.center();

            // get solution (eq 24)
            const auto solution = tgv.u( center_position * 2 * M_PI, finalTime);

            // compute error vectors (eq 25)
            const double err_ux = u[cellIdx, 0] - solution[0];
            const double err_uy = u[cellIdx, 1] - solution[1];
            const double err_uz = u[cellIdx, 2] - solution[2];
            // store errors
            velocityError[cellIdx, 0] = err_ux;
            velocityError[cellIdx, 1] = err_uy;
            velocityError[cellIdx, 2] = err_uz;
            // accumulate L2 norms (eq 26)
            error_L2_ux += err_ux * err_ux;
            error_L2_uy += err_uy * err_uy;
            norm_L2_ux += solution[0] * solution[0];
            norm_L2_uy += solution[1] * solution[1];
        }

        // compute final L2 errors (eq 26)
        const double relative_error_ux = std::sqrt(error_L2_ux / norm_L2_ux);
        const double relative_error_uy = std::sqrt(error_L2_uy / norm_L2_uy);

        std::println(std::clog, "Validation Results:");
        std::println(std::clog, "  Relative L2 Error in u_x: {:.6e}", relative_error_ux);
        std::println(std::clog, "  Relative L2 Error in u_y: {:.6e}", relative_error_uy);

        // write errors to file
        const std::filesystem::path errorOutputPath = outputFile / "errors.txt";
        std::ofstream errorFile(errorOutputPath);
        if (!errorFile) {
            throw std::runtime_error("Failed to open error output file");
        }
        errorFile << std::format("{:.16e} {:.16e}\n", relative_error_ux, relative_error_uy);
        errorFile.close();
        std::println(std::clog, "Errors written to {}", errorOutputPath.string());

        // export final results
        const std::filesystem::path outputPath = outputFile / "tgv_results.vtkhdf";
        auto exporter = oktal::io::vtk::exportCellGrid(grid,outputPath);
        exporter.writeGridVector<double, 1>("density", rho.const_view());
        exporter.writeGridVector<double, 3>("velocity", u.const_view());
        exporter.writeGridVector<double, 3>("velocity_error", velocityError.const_view());

        std::println(std::clog, "Final Results exported to {}", outputPath.string());

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