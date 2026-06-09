#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <map>
#include <utility>

#include "coot-utils/atom-vertex.hh"
#include "coot-utils/map-index.hh"
#include "geometry/protein-geometry.hh"

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__atom_tree_t__construct_internal/gemmi/function.hh"

namespace coot {

class atom_tree_t {
public:
    std::vector<atom_vertex> atom_vertex_vec;
    std::map<std::string, map_index_t> name_to_index;
    int made_from_minimol_residue_flag;
    int n_selected_atoms;
    const void *atom_selection;

    void atom_tree_t_gemmi(
        const dictionary_residue_restraints_t &rest,
        gemmi::Residue *res,
        const std::string &altconf) {
        made_from_minimol_residue_flag = 0;
        n_selected_atoms = 0;
        atom_selection = nullptr;
        std::vector<std::pair<int, int>> bonds;
        construct_internal_gemmi(rest, res, altconf, bonds, name_to_index, atom_vertex_vec);
    }
};

} // namespace coot
