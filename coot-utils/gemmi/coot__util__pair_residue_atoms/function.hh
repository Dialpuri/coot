#pragma once

#include <vector>
#include <utility>
#include <string>
#include <gemmi/model.hpp>

namespace coot {
namespace util {

std::vector<std::pair<int, int>>
pair_residue_atoms_gemmi(const gemmi::Residue *a_residue,
                         const gemmi::Residue *b_residue) {

   std::vector<std::pair<int, int>> pv;

   const auto &atoms1 = a_residue->atoms;
   const auto &atoms2 = b_residue->atoms;

   for (int i = 0; i < (int)atoms1.size(); i++) {
      std::string atn1 = atoms1[i].name;
      std::string alt1(1, atoms1[i].altloc);
      for (int j = 0; j < (int)atoms2.size(); j++) {
         std::string atn2 = atoms2[j].name;
         std::string alt2(1, atoms2[j].altloc);
         if (atn1 == atn2) {
            if (alt1 == alt2) {
               std::pair<int, int> p(i, j);
               pv.push_back(p);
               break;
            }
         }
      }
   }
   return pv;
}

} // namespace util
} // namespace coot