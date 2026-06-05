#pragma once

#include <gemmi/model.hpp>
#include <gemmi/math.hpp>
#include <vector>
#include <string>
#include "coot-utils/map-index.hh"
#include "coot-utils/atom-vertex.hh"

namespace coot {

double atom_tree_t_rotate_about_gemmi(
    const std::vector<atom_vertex>& atom_vertex_vec,
    const std::vector<const gemmi::Atom*>& atom_selection,
    gemmi::Residue* residue,
    int index2,
    int index3,
    double angle,
    bool reversed_flag);

} // namespace coot