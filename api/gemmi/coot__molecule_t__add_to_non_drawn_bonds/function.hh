#pragma once
#include <string>
#include <set>
#include <gemmi/model.hpp>

namespace coot {

// Port of molecule_t::add_to_non_drawn_bonds to gemmi
// Since gemmi doesn't have UDD handles or the same atom indexing system,
// this version takes the structure and a mutable reference to the set of indices to modify
// The indices are relative to the flattened atom list in the model
void molecule_t_add_to_non_drawn_bonds_gemmi(const gemmi::Structure& st, 
                                             std::set<int>& no_bonds_to_these_atom_indices,
                                             const std::string& atom_selection_cid);

}