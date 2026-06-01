#pragma once

#include <string>
#include <vector>
#include <utility>
#include <gemmi/model.hpp>
#include "geometry/residue-and-atom-specs.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_residue__ffd79d/gemmi/function.hh"
#include "clipper/core/coords.h"

namespace coot {
namespace util {

inline std::pair<clipper::Coord_orth, clipper::Coord_orth>
extents_gemmi(gemmi::Structure &st,
              const std::vector<residue_spec_t> &specs) {

   float most_x = -99999;
   float most_y = -99999;
   float most_z = -99999;
   float least_x = 99999;
   float least_y = 99999;
   float least_z = 99999;
   for (unsigned int ispec=0; ispec<specs.size(); ispec++) {
      gemmi::Residue *residue_p = get_residue_gemmi(specs[ispec], st);
      if (residue_p) {
         for (const gemmi::Atom &atom : residue_p->atoms) {
            if (atom.pos.x < least_x) least_x = static_cast<float>(atom.pos.x);
            if (atom.pos.y < least_y) least_y = static_cast<float>(atom.pos.y);
            if (atom.pos.z < least_z) least_z = static_cast<float>(atom.pos.z);
            if (atom.pos.x >  most_x)  most_x  = static_cast<float>(atom.pos.x);
            if (atom.pos.y >  most_y)  most_y  = static_cast<float>(atom.pos.y);
            if (atom.pos.z >  most_z)  most_z  = static_cast<float>(atom.pos.z);
         }
      }
   }
   clipper::Coord_orth p1( most_x,  most_y,  most_z);
   clipper::Coord_orth p2(least_x, least_y, least_z);
   return std::pair<clipper::Coord_orth, clipper::Coord_orth> (p2, p1);
}

} // namespace util
} // namespace coot