#pragma once
#include <string>
#include <vector>
#include <gemmi/model.hpp>

namespace coot {
namespace molecule_t {

// Adjust occupancies of other atoms in the same residue with same name but different altLoc
void adjust_occupancy_other_residue_atoms_gemmi(
    gemmi::Atom* atom,
    gemmi::Residue* residue,
    bool force_sum_1_flag);

}
}