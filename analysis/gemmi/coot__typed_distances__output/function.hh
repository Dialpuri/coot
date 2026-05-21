#pragma once
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <gemmi/model.hpp>

namespace coot {
namespace typed_distances {

// Helper function to convert gemmi residue to a string representation
inline std::string residue_to_string(const gemmi::Residue& res) {
    return res.name + " " + std::to_string(res.seqid.num.value);
}

// Port of coot::typed_distances::output() to gemmi
// This function outputs a windowed residue table showing distance statistics
void output_gemmi(
    const std::map<const gemmi::Residue*, std::map<int, std::vector<float>>>& residue_distances_map,
    const std::map<const gemmi::Residue*, std::vector<const gemmi::Residue*>>& residues_within_window
);

} // namespace typed_distances
} // namespace coot