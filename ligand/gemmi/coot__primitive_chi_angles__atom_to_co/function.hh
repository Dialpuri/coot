#pragma once
#include <clipper/core/coords.h>
#include <gemmi/model.hpp>

namespace coot {
namespace primitive_chi_angles {

// Ported from: clipper::Coord_orth atom_to_co(mmdb::Atom *at) const
// Original: returns clipper::Coord_orth with atom coordinates
clipper::Coord_orth atom_to_co_gemmi(const gemmi::Atom& atom);

} // namespace primitive_chi_angles
} // namespace coot