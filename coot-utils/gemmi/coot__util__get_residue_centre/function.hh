#pragma once
#include <utility>
#include <clipper/core/coords.h>
#include <gemmi/model.hpp>

namespace coot { namespace util {

// Port of coot::util::get_residue_centre from MMDB to gemmi
// Returns (status, centre) where status is true if residue has atoms
inline std::pair<bool, clipper::Coord_orth>
get_residue_centre_gemmi(const gemmi::Residue& residue) {
   bool status = false;
   clipper::Coord_orth centre(0,0,0);

   int n_residue_atoms = residue.atoms.size();
   if (n_residue_atoms > 0) {
      status = true;
      for (const auto& atom : residue.atoms) {
         clipper::Coord_orth pt(atom.pos.x, atom.pos.y, atom.pos.z);
         centre += pt;
      }
      double scale = 1.0 / static_cast<double>(n_residue_atoms);
      centre = scale * centre;
   }
   return std::make_pair(status, centre);
}

}}