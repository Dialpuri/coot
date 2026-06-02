#pragma once

#include <vector>
#include <utility>
#include <gemmi/model.hpp>

namespace coot {

// gemmi port of restraints_container_t::make_helix_pseudo_bond_restraints_from_res_vec_auto
// Returns the number of helical pseudo-bond restraints that would be added.
int make_helix_pseudo_bond_restraints_from_res_vec_auto_gemmi(
    const std::vector<std::pair<bool, gemmi::Residue*>>& residues_vec);

} // namespace coot