#pragma once

#include <string>
#include <gemmi/model.hpp>

namespace coot {

class bonded_pair_t {
public:
   void delete_atom_gemmi(gemmi::Residue *res, const std::string &atom_name);
};

} // namespace coot

inline void coot::bonded_pair_t::delete_atom_gemmi(gemmi::Residue *res, const std::string &atom_name) {
   for (auto it = res->atoms.begin(); it != res->atoms.end(); ++it) {
      if (it->name == atom_name) {
         res->atoms.erase(it);
         return;
      }
   }
}