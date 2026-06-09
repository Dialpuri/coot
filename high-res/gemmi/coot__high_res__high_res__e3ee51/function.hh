#pragma once

#include "mini-mol/mini-mol.hh"
#include "clipper/core/coords.h"
#include <gemmi/model.hpp>

namespace coot {
namespace high_res {

// Gemmi port of coot::high_res::high_res (the constructor).
// Takes a molecule, a centre, and a gemmi::Structure (for neighbor search),
// returns the globular molecule.
coot::minimol::molecule high_res_gemmi(
    const coot::minimol::molecule &mol,
    const clipper::Coord_orth &given_centre,
    const gemmi::Structure &structure);

} // namespace high_res
} // namespace coot
