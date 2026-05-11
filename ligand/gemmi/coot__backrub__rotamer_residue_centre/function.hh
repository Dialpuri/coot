#pragma once

#include <clipper/clipper.h>
#include <gemmi/model.hpp>

namespace coot {
namespace backrub {

// Ported from coot::backrub::rotamer_residue_centre()
// Calculates the centroid (mean position) of all atoms in a residue
clipper::Coord_orth rotamer_residue_centre_gemmi(const gemmi::Residue& residue);

} // namespace backrub
} // namespace coot