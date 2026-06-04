#pragma once

#include <clipper/core/coords.h>
#include <gemmi/model.hpp>
#include <utility>
#include <vector>

namespace coot {
namespace util {

// Existing coot utility — not MMDB, so we call it directly
clipper::RTop_orth make_rtop_orth_for_jiggle_atoms(
    float jiggle_trans_scale_factor, float annealing_factor);

std::pair<clipper::RTop_orth, std::vector<gemmi::Atom>>
jiggle_atoms_gemmi(const std::vector<gemmi::Atom>& atoms,
                   const clipper::Coord_orth& centre_pt,
                   float jiggle_trans_scale_factor,
                   float annealing_factor) {

  std::vector<gemmi::Atom> new_atoms(atoms.size());
  clipper::RTop_orth rtop = make_rtop_orth_for_jiggle_atoms(jiggle_trans_scale_factor, annealing_factor);
  for (unsigned int i = 0; i < atoms.size(); i++) {
    clipper::Coord_orth pt_rel(atoms[i].pos.x - centre_pt.x(),
                               atoms[i].pos.y - centre_pt.y(),
                               atoms[i].pos.z - centre_pt.z());
    clipper::Coord_orth new_pt = pt_rel.transform(rtop);
    new_pt += centre_pt;
    new_atoms[i] = atoms[i];
    new_atoms[i].pos.x = new_pt.x();
    new_atoms[i].pos.y = new_pt.y();
    new_atoms[i].pos.z = new_pt.z();
  }
  return std::pair<clipper::RTop_orth, std::vector<gemmi::Atom>>(rtop, new_atoms);
}

} // namespace util
} // namespace coot