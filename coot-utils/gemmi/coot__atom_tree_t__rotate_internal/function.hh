#pragma once

#include <gemmi/model.hpp>
#include <gemmi/math.hpp>
#include <vector>
#include "coot-utils/map-index.hh"

namespace coot {

class atom_tree_t {
public:
  void rotate_internal_gemmi(
    std::vector<gemmi::Atom> &atoms,
    std::vector<coot::map_index_t> moving_atom_indices,
    const gemmi::Vec3 &dir,
    const gemmi::Vec3 &base_atom_pos,
    double angle);
};

} // namespace coot