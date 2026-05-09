#pragma once
#include <vector>
#include <string>
#include <gemmi/model.hpp>
#include "coot/api/moved-atom.hh"

namespace coot {
namespace molecule_t {

// Port of new_positions_for_residue_atoms from MMDB to gemmi
// Returns number of atoms processed, or -1 on error
int new_positions_for_residue_atoms_gemmi(
    gemmi::Model& model,
    const std::string& residue_cid,
    const std::vector<coot::api::moved_atom_t>& moved_atoms);

} // namespace molecule_t
} // namespace coot