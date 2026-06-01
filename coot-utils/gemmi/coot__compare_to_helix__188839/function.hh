#pragma once

#include <vector>
#include <string>
#include <cmath>

#include <clipper/core/coords.h>

#include <gemmi/model.hpp>
#include "coot-utils/helix-like.hh"

namespace coot {

helical_results_t compare_to_helix_gemmi(
    const std::vector<gemmi::Residue *> &helical_residues,
    const std::vector<clipper::Coord_orth> &alpha_ref_positions) {

  double sum_delta_lim = 3.6;

  helical_results_t r;
  std::vector<clipper::Coord_orth> match_set(16); // 4 x 4

  if (helical_residues.size() == 4) {
    int n_found = 0;
    for (unsigned int i = 0; i < 4; i++) {
      gemmi::Residue *residue_p = helical_residues[i];
      for (const gemmi::Atom &at : residue_p->atoms) {
        int idx = -1;
        std::string atom_name = at.name;
        if (atom_name == "N") idx = 0;
        if (atom_name == "CA") idx = 1;
        if (atom_name == "C") idx = 2;
        if (atom_name == "O") idx = 3;
        if (idx != -1) {
          int idx_match_set = static_cast<int>(i) * 4 + idx;
          clipper::Coord_orth co(at.pos.x, at.pos.y, at.pos.z);
          match_set[idx_match_set] = co;
          n_found++;
        }
      }
    }

    if (n_found == 16) {
      clipper::RTop_orth rtop(alpha_ref_positions, match_set);
      double sum_delta = 0.0;
      for (unsigned int ii = 0; ii < 16; ii++) {
        clipper::Coord_orth moved_pos = rtop * alpha_ref_positions[ii];
        double dd = clipper::Coord_orth(match_set[ii] - moved_pos).lengthsq();
        sum_delta += std::sqrt(dd);
      }
      if (sum_delta < sum_delta_lim) {
        r.is_alpha_helix_like = true;
      }
      r.sum_delta = sum_delta;
    }
  }
  return r;
}

helical_results_t compare_to_helix_gemmi(
    const std::vector<gemmi::Residue *> &helical_residues) {

  helical_results_t r;
  std::vector<clipper::Coord_orth> ref_pos = coot::alpha_helical_reference_positions();

  if (helical_residues.size() == 4)
    r = compare_to_helix_gemmi(helical_residues, ref_pos);
  return r;
}

} // namespace coot