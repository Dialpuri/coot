#pragma once

#include <clipper/clipper.h>
#include <gemmi/model.hpp>
#include <vector>

namespace coot {
namespace util {
namespace fffear_search {

int fill_nxmap_gemmi(clipper::NXmap<float>& nxmap,
                     const std::vector<const gemmi::Atom*>& atoms,
                     const clipper::Coord_orth& mid_point);

} // namespace fffear_search
} // namespace util
} // namespace coot