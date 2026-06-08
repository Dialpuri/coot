#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <map>
#include <utility>

#include "geometry/protein-geometry.hh"
#include "coot-utils/atom-vertex.hh"
#include "coot-utils/map-index.hh"

// Forward declarations for the _gemmi port functions
namespace coot {

std::map<std::string, map_index_t> fill_name_map_gemmi(
    const gemmi::Residue &res,
    const std::string &altconf
);

bool fill_atom_vertex_vec_gemmi(
    const dictionary_residue_restraints_t &rest,
    gemmi::Residue *res,
    const std::string &altconf,
    std::map<std::string, map_index_t> &name_to_index,
    std::vector<atom_vertex> &atom_vertex_vec
);

void construct_internal_gemmi(
    const dictionary_residue_restraints_t &rest,
    gemmi::Residue *res,
    const std::string &altconf,
    std::vector<std::pair<int, int>> &bonds_out,
    std::map<std::string, map_index_t> &name_to_index_out,
    std::vector<atom_vertex> &atom_vertex_vec_out
);

} // namespace coot
