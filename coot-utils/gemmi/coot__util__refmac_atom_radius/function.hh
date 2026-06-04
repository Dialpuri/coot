#pragma once

#include <gemmi/model.hpp>
#include <cmath>

namespace coot {
namespace util {

inline double refmac_atom_radius_gemmi(gemmi::Atom const& at) {
    double u = at.b_iso / (8 * M_PI * M_PI);
    double v = 2 * u;
    v = std::max(v, 0.38); // don't have tiny atoms
    v += 0.6;
    v *= 4;
    return std::sqrt(v);
}

} // namespace util
} // namespace coot