#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <iostream>

namespace coot {

/// Port of coot::read_amber_trajectory to gemmi.
/// Like the original, this stub always returns nullptr because
/// NetCDF support is not compiled.
///
/// @param topology_mol      Pointer to the topology structure (may be nullptr).
/// @param trajectory_file_name  Path to the AMBER NetCDF trajectory file.
/// @param start_frame       First frame to read (1-based).
/// @param end_frame         Last frame to read (1-based).
/// @param stride            Frame stride (1 = every frame).
/// @return nullptr (NetCDF not available).
inline gemmi::Structure* read_amber_trajectory_gemmi(
    gemmi::Structure* topology_mol,
    const std::string& trajectory_file_name,
    int start_frame,
    int end_frame,
    int stride) {
    (void)topology_mol;
    (void)trajectory_file_name;
    (void)start_frame;
    (void)end_frame;
    (void)stride;
    std::cout << "WARNING:: read_amber_trajectory: NetCDF support not compiled" << std::endl;
    return nullptr;
}

}  // namespace coot