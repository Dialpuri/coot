#pragma once
#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include <cmath>

namespace coot {
namespace reduce {

std::pair<clipper::Coord_orth, clipper::Coord_orth>
position_pair_by_bisection_gemmi(const gemmi::Atom* at_1,
                                  const gemmi::Atom* at_2,
                                  const gemmi::Atom* at_3,
                                  double bl,
                                  double alpha);

} // namespace reduce
} // namespace coot