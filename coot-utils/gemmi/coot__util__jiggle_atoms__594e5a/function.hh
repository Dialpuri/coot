#pragma once

#include <gemmi/model.hpp>
#include <utility>
#include <vector>

#include "clipper/core/coords.h"

namespace coot {
namespace util {

std::pair<clipper::RTop_orth, std::vector<gemmi::Atom>>
jiggle_atoms_gemmi(
    std::vector<gemmi::Atom const *> const &atoms,
    clipper::Coord_orth const &centre_pt,
    float jiggle_trans_scale_factor,
    float annealing_factor);

} // namespace util
} // namespace coot