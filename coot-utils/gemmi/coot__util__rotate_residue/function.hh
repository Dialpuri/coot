#pragma once

#include <gemmi/model.hpp>
#include <clipper/core/coords.h>

// Include the definition of rotate_around_vector
#include "coot-utils/coot-coord-utils.hh"

namespace coot {
namespace util {

// Rotate all atoms in a residue around an axis defined by `direction`,
// passing through `origin_shift`, by `angle` radians.
//
// Ported from MMDB version: iterates residue.atoms instead of using
// GetAtomTable, and skips the isTer() check which has no gemmi equivalent
// (TER records don't produce atoms in gemmi's data model).
inline void rotate_residue_gemmi(
    gemmi::Residue* residue_p,
    const clipper::Coord_orth& direction,
    const clipper::Coord_orth& origin_shift,
    double angle) {

   if (residue_p) {
      for (gemmi::Atom& at : residue_p->atoms) {
         clipper::Coord_orth pt(at.pos.x, at.pos.y, at.pos.z);
         clipper::Coord_orth pt_new = rotate_around_vector(direction, pt, origin_shift, angle);
         at.pos.x = pt_new.x();
         at.pos.y = pt_new.y();
         at.pos.z = pt_new.z();
      }
   }
}

} // namespace util
} // namespace coot