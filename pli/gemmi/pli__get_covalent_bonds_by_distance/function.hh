#pragma once

#include <gemmi/model.hpp>
#include "pli/protein-ligand-interactions.hh"
#include "geometry/residue-and-atom-specs.hh"
#include <string>
#include <vector>

namespace pli {

std::vector<pli::fle_ligand_bond_t>
get_covalent_bonds_by_distance_gemmi(const gemmi::Model& model,
                                     const coot::residue_spec_t& ligand_spec);

} // namespace pli