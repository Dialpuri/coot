#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <string>
#include <map>
#include "coot-utils/atom-vertex.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__atom_tree_t__fill_atom_vertex_vec_using_contacts/gemmi/function.hh"

namespace coot {

void atom_tree_t_gemmi(
    const std::vector<std::vector<int>>& contact_indices,
    int base_atom_index,
    gemmi::Residue* res,
    const std::string& altconf,
    std::map<std::string, int>& name_to_index,
    std::vector<atom_vertex>& atom_vertex_vec);

} // namespace coot