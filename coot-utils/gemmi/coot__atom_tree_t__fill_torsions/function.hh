#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>

#include "coot-utils/atom-tree.hh"
#include "coot-utils/atom-vertex.hh"
#include "geometry/protein-geometry.hh"

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__atom_tree_t__get_atom_index_quad/gemmi/function.hh"

namespace coot {

// gemmi port of atom_tree_t::fill_torsions
// Follows the same pattern as construct_internal_gemmi: the atom_vertex_vec
// is passed by reference (it is protected in atom_tree_t).
bool fill_torsions_gemmi(
    std::vector<coot::atom_vertex> &atom_vertex_vec,
    const dictionary_residue_restraints_t &rest,
    gemmi::Residue *res,
    const std::string &altconf);

} // namespace coot
