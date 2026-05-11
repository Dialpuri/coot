#pragma once
#include <gemmi/model.hpp>
#include <clipper/core/coords.h>

namespace coot {

// Port of coot::co from MMDB to gemmi
// Takes a gemmi::Atom pointer and returns clipper::Coord_orth with its coordinates
inline clipper::Coord_orth co_gemmi(const gemmi::Atom* at) {
    return clipper::Coord_orth(at->pos.x, at->pos.y, at->pos.z);
}

} // namespace coot