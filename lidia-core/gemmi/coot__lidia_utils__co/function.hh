#pragma once

#include <clipper/core/coords.h>
#include <gemmi/model.hpp>

namespace coot {
namespace lidia_utils {

inline clipper::Coord_orth
co_gemmi(const gemmi::Atom *at) {
    return clipper::Coord_orth(at->pos.x, at->pos.y, at->pos.z);
}

} // namespace lidia_utils
} // namespace coot