#pragma once

#include <clipper/core/xmap.h>
#include <clipper/core/cell.h>
#include <clipper/core/spacegroup.h>
#include <clipper/core/coords.h>
#include <gemmi/model.hpp>
#include <vector>

namespace coot {
namespace util {

clipper::Xmap<float> make_map_mask_gemmi(
    const clipper::Spacegroup &space_group,
    const clipper::Cell &cell,
    const clipper::Grid_sampling &grid_sampling,
    std::vector<const gemmi::Atom*> atoms,
    float radius,
    float smooth);

} // namespace util
} // namespace coot