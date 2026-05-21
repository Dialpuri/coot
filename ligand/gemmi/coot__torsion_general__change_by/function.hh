#pragma once
#include <vector>
#include <gemmi/model.hpp>

namespace coot {
namespace torsion_general {

// Port of change_by function - rotates atoms around a bond by given torsion angle (in degrees)
// Returns 0 on success, non-zero on failure
int change_by_gemmi(double diff_degrees,
                    std::vector<gemmi::Position>& coords,
                    int atom1_idx, int atom2_idx,
                    const std::vector<int>& children_of_atom2);

} // namespace torsion_general
} // namespace coot