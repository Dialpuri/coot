#pragma once
#include <clipper/core/coords.h>
#include <gemmi/model.hpp>

namespace coot {

// Ported from: clipper::Coord_orth monomer_utils::atom_to_co(mmdb::Atom *at) const
// Original: return clipper::Coord_orth(at->x, at->y, at->z);
inline clipper::Coord_orth atom_to_co_gemmi(const gemmi::Atom& atom) {
    return clipper::Coord_orth(atom.pos.x, atom.pos.y, atom.pos.z);
}

} // namespace coot