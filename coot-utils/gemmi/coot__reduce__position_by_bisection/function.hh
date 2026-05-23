#pragma once
#include <clipper/core/coords.h>
#include <gemmi/model.hpp>

namespace coot {
namespace reduce {

clipper::Coord_orth position_by_bisection_gemmi(const gemmi::Atom* at_1,
                                                const gemmi::Atom* at_2,
                                                const gemmi::Atom* at_3,
                                                double bl);

} // namespace reduce
} // namespace coot