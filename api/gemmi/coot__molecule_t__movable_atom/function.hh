#pragma once

#include <gemmi/model.hpp>
#include <iostream>

namespace coot {

inline bool movable_atom_gemmi(const gemmi::Atom *atom, bool replace_coords_with_zero_occ_flag) {
   if (!atom) {
      std::cout << "ERROR:: null atom in movable_atom()" << std::endl;
      return false;
   }

   bool m = true;

   if ((atom->occ < 0.0001) && (atom->occ > -0.0001))
      if (!replace_coords_with_zero_occ_flag)
         m = false;
   return m;
}

} // namespace coot