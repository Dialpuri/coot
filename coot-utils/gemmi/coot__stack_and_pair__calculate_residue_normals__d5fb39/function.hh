#pragma once

#include <map>
#include <set>
#include <utility>
#include <vector>
#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__stack_and_pair__get_base_normal/gemmi/function.hh"

namespace coot {
namespace stack_and_pair {

std::map<const gemmi::Residue*, clipper::Coord_orth>
calculate_residue_normals_gemmi(std::vector<const gemmi::Residue*> residues_vec) {

   std::set<const gemmi::Residue*> done_res;
   std::map<const gemmi::Residue*, clipper::Coord_orth> m;
   for (std::size_t i = 0; i < residues_vec.size(); i++) {
      const gemmi::Residue* r = residues_vec[i];
      if (r == nullptr) continue;
      if (done_res.find(r) == done_res.end()) {
         std::pair<bool, clipper::Coord_orth> bn = get_base_normal_gemmi(r);
         if (bn.first) {
            m[r] = bn.second;
         }
         done_res.insert(r);
      }
   }
   return m;
}

} // namespace stack_and_pair
} // namespace coot