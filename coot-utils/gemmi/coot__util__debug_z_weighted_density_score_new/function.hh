#pragma once

#include <gemmi/model.hpp>
#include <clipper/core/xmap.h>
#include <vector>
#include <utility>
#include <iostream>
#include "coot-utils/coot-map-utils.hh"

namespace coot {
namespace util {

void debug_z_weighted_density_score_new_gemmi(
    const std::vector<std::pair<gemmi::Atom *, float>> &atom_atom_number_pairs,
    const clipper::Xmap<float> &map) {

   float sum_d = 0;
   for (unsigned int iat = 0; iat < atom_atom_number_pairs.size(); iat++) {
      gemmi::Atom *at = atom_atom_number_pairs[iat].first;
      if (!at) continue;
      clipper::Coord_orth co(at->pos.x, at->pos.y, at->pos.z);
      float d = coot::util::density_at_point(map, co);
      float w = atom_atom_number_pairs[iat].second;
      sum_d += d * w;
      std::cout << "debug score " << iat << " atom_spec "
                << at->name
                << " pos " << at->pos.x << " " << at->pos.y << " " << at->pos.z
                << " weight: " << w << " density:" << d << " running sum " << sum_d << std::endl;
   }
   std::cout << "debug:: debug_z_weighted_density_score_new(): total: " << sum_d << std::endl;
}

} // namespace util
} // namespace coot
