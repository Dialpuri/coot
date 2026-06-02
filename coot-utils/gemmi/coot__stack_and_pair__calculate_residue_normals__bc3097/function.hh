#pragma once

#include <map>
#include <utility>
#include <vector>
#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__stack_and_pair__get_base_normal/gemmi/function.hh"

namespace coot {
namespace stack_and_pair {

std::map<const gemmi::Residue*, clipper::Coord_orth>
calculate_residue_normals_gemmi(const std::vector<std::pair<bool, const gemmi::Residue*>>& residues_vec) {

   std::map<const gemmi::Residue*, clipper::Coord_orth> m;
   for (std::size_t i = 0; i < residues_vec.size(); i++) {
      const gemmi::Residue* r = residues_vec[i].second;
      if (r) {
         std::pair<bool, clipper::Coord_orth> normal = get_base_normal_gemmi(r);
         if (normal.first) {
            m[r] = normal.second;
         }
      }
   }
   return m;
}

} // namespace stack_and_pair
} // namespace coot