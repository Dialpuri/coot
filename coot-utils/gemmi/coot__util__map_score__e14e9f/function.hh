#pragma once

#include <vector>
#include <gemmi/model.hpp>
#include "clipper/core/xmap.h"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"
#include "coot-utils/coot-map-utils.hh"

namespace coot {
namespace util {

float map_score_gemmi(const std::vector<const gemmi::Atom*> &atom_selection,
                      const clipper::Xmap<float> &xmap,
                      short int with_atomic_weighting) {
    float f = 0.0;
    for (const gemmi::Atom* at : atom_selection) {
        clipper::Coord_orth co = coot::co_gemmi(at);
        float f1 = density_at_point(xmap, co);
        f1 *= at->occ;
        f += f1;
    }
    return f;
}

} // namespace util
} // namespace coot