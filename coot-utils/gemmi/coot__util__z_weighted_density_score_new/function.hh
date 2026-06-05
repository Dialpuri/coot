#pragma once

#include <vector>
#include <utility>
#include <clipper/core/coords.h>
#include <clipper/core/xmap.h>
#include <gemmi/model.hpp>

#include "coot-utils/coot-map-utils.hh"

namespace coot {
namespace util {

float z_weighted_density_score_new_gemmi(
    const std::vector<std::pair<gemmi::Atom*, float>>& atom_weight_pairs,
    const clipper::Xmap<float>& map
) {
    float sum_d = 0;
    for (const auto& p : atom_weight_pairs) {
        gemmi::Atom* at = p.first;
        clipper::Coord_orth co(at->pos.x, at->pos.y, at->pos.z);
        float d = coot::util::density_at_point(map, co) * p.second;
        sum_d += d;
    }
    return sum_d;
}

} // namespace util
} // namespace coot