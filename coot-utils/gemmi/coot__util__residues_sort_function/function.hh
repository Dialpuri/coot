#pragma once

#include <gemmi/model.hpp>
#include <string>

namespace coot {
namespace util {

bool residues_sort_function_gemmi(gemmi::CRA cra1, gemmi::CRA cra2) {

   if (!cra1.residue) {
      if (!cra2.residue) {
         return false;
      } else {
         return true;
      }
   } else {
      if (!cra2.residue) {
         return false;
      }
   }

   if (cra1.chain < cra2.chain) {
      return true;
   } else {
      if (cra1.chain > cra2.chain) {
         return false;
      } else {
         if (cra1.residue->seqid.num.value < cra2.residue->seqid.num.value) {
            return true;
         } else {
            if (cra1.residue->seqid.num.value > cra2.residue->seqid.num.value) {
               return false;
            } else {
               auto norm = [](char c) -> std::string {
                  return (c == '\0') ? "" : std::string(1, c);
               };
               std::string inscode_1 = norm(cra1.residue->seqid.icode);
               std::string inscode_2 = norm(cra2.residue->seqid.icode);
               return (inscode_1 < inscode_2);
            }
         }
      }
   }
}

} // namespace util
} // namespace coot