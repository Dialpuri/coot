#pragma once

#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include <vector>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot {

class atom_overlaps_container_t {

 public:
  bool is_inside_another_atom_to_which_its_bonded_gemmi(
      int atom_idx,
      const gemmi::Atom *at,
      const clipper::Coord_orth &pt_at_surface,
      const std::vector<int> &bonded_neighb_indices,
      const gemmi::Atom **atom_selection,
      const std::vector<double> &neighb_atom_radius) {

    bool r = false;
    double r_1 = neighb_atom_radius[atom_idx];

    for (unsigned int i = 0; i < bonded_neighb_indices.size(); i++) {
      const gemmi::Atom *clash_neighb = atom_selection[bonded_neighb_indices[i]];
      clipper::Coord_orth pt_clash_neigh = coot::co_gemmi(clash_neighb);
      double r_2 = neighb_atom_radius[bonded_neighb_indices[i]];
      double r_2_sqrd = r_2 * r_2;
      double d_sqrd = (pt_at_surface - pt_clash_neigh).lengthsq();
      if (d_sqrd < r_2_sqrd) {
        r = true;
        break;
      }
    }

    return r;
  }
};

}  // namespace coot