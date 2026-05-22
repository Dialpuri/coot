#pragma once
#include <vector>
#include <gemmi/model.hpp>
#include "coot-utils/atom-vertex.hh"

namespace coot {

// Gemmi port of MMDB's fill_atom_vertex_vec_using_contacts
// Ported from:
//   bool fill_atom_vertex_vec_using_contacts(const std::vector<std::vector<int> > &contact_indices, int base_atom_index);
//   which used mmdb::Residue* and called GetAtomTable()
//
// Now uses gemmi::Residue and extracts atoms directly from it.
// Calls the gemmi port of fill_atom_vertex_vec_using_contacts_by_atom_selection.
bool fill_atom_vertex_vec_using_contacts_gemmi(const std::vector<std::vector<int> > &contact_indices,
                                               int base_atom_index,
                                               const std::vector<gemmi::Atom>& residue_atoms,
                                               std::vector<atom_vertex>& atom_vertex_vec);

} // namespace coot