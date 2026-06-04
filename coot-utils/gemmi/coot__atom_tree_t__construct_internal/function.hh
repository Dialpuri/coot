#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <map>
#include <utility>

#include "coot-utils/atom-tree.hh"
#include "coot-utils/map-index.hh"
#include "geometry/protein-geometry.hh"

namespace coot {

void construct_internal_gemmi(
    const dictionary_residue_restraints_t& rest,
    gemmi::Residue* res,
    const std::string& altconf,
    std::vector<std::pair<int,int>>& bonds_out,
    std::map<std::string, map_index_t>& name_to_index_out,
    std::vector<atom_vertex>& atom_vertex_vec_out);

} // namespace coot