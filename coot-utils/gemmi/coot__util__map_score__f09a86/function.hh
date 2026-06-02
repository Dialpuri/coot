#pragma once

#include <vector>
#include "clipper/core/xmap.h"
#include <gemmi/model.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"
#include "coot-utils/coot-map-utils.hh"

namespace coot {
namespace util {

inline float map_score_gemmi(std::vector<const gemmi::Atom *> atoms,
                             const clipper::Xmap<float> &xmap) {
  float f = 0.0f;
  for (unsigned int i = 0; i < atoms.size(); i++) {
    if (atoms[i]) {
      float f1 = density_at_point(xmap, co_gemmi(atoms[i]));
      f1 *= atoms[i]->occ;
      f += f1;
    }
  }
  return f;
}

} // namespace util
} // namespace coot