#pragma once

#include <gemmi/model.hpp>
#include <map>
#include <vector>
#include <string>
#include "coot-utils/atom-vertex.hh"
#include "coot-utils/map-index.hh"
#include "geometry/protein-geometry.hh"

namespace coot {

// gemmi port of coot::atom_tree_t::fill_atom_vertex_vec
// Freestanding because atom_tree_t's data members are protected/private.
bool fill_atom_vertex_vec_gemmi(
    const dictionary_residue_restraints_t& rest,
    gemmi::Residue* res,
    const std::string& altconf,
    std::map<std::string, map_index_t>& name_to_index,
    std::vector<atom_vertex>& atom_vertex_vec);

// Helper: build name_to_index and atom_vertex_vec from a gemmi::Residue
// (mirrors the constructor logic from the existing atom_tree_t_gemmi port)
void atom_tree_t_gemmi(
    const std::vector<std::vector<int> >& contact_indices,
    int base_atom_index,
    gemmi::Residue* res,
    const std::string& altconf,
    std::map<std::string, int>& name_to_index_out,
    std::vector<atom_vertex>& atom_vertex_vec_out);

} // namespace coot