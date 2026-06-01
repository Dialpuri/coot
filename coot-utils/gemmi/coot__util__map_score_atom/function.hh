#pragma once

#include <gemmi/model.hpp>
#include <clipper/core/xmap.h>

namespace coot {
namespace util {

// gemmi port of coot::util::map_score_atom
// Returns the electron density at the atom's position from the given map.
// Returns 0.0 for a null atom pointer.
float map_score_atom_gemmi(const gemmi::Atom *atom,
                           const clipper::Xmap<float> &xmap);

} // namespace util
} // namespace coot