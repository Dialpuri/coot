#pragma once
#include <string>
#include <gemmi/model.hpp>

namespace coot {
namespace molecule_t {

// Port of coot::molecule_t::set_residue_to_rotamer_move_atoms to gemmi
// Copies atom coordinates from moving_res to res for atoms with matching names and altLocs
// Returns 1 if any atoms were moved, 0 otherwise
int set_residue_to_rotamer_move_atoms_gemmi(gemmi::Residue& res, gemmi::Residue& moving_res);

} // namespace molecule_t
} // namespace coot