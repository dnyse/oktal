#include <iostream>
#include <filesystem>
#include <string>
#include <format>
#include <print>
#include <span>
#include "oktal/io/VtkExport.hpp"
#include "oktal/octree/CellOctree.hpp"


int main(int argc, char* argv[]) {

    try {

        const std::span<char*> args(argv, static_cast<std::size_t>(argc));

        if (argc != 3) {
            std::println(std::clog, "Usage: {} <filepath> <descriptor>", args[0]); 
            return 1;
        }

        const std::filesystem::path filepath = args[1];
        const std::string descriptor = args[2];    

        const oktal::CellOctree octree = oktal::CellOctree::fromDescriptor(descriptor);
        
        oktal::io::vtk::exportOctree(octree, filepath);

        std::println("Successfully created {}", filepath.string());

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