#pragma once

#include <utility>
#include <string>
#include <gemmi/model.hpp>
#include <clipper/core/coords.h>

namespace coot { namespace util {

std::pair<bool, clipper::Coord_orth>
get_CA_position_in_residue_gemmi(const gemmi::Residue &res) {
   for (const auto &at : res.atoms) {
      if (at.name == "CA")
         return std::make_pair(true, clipper::Coord_orth(at.pos.x, at.pos.y, at.pos.z));
   }
   return std::make_pair(false, clipper::Coord_orth(0,0,0));
}

}} // namespace coot::util