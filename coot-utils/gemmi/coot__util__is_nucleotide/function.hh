#pragma once
#include <string>
#include <gemmi/model.hpp>

namespace coot { namespace util {

// Return 1 if the residue is a nucleotide (RNA or DNA), 0 otherwise
// Matches the original MMDB function logic
inline short int
is_nucleotide_gemmi(const gemmi::Residue &r) {
   std::string rn = r.name;
   if (rn == "A"  || rn == "C"  || rn == "G"  || rn == "U"  || rn == "T"  ||
       rn == "DA" || rn == "DC" || rn == "DG" || rn == "DT" ||
       rn == "Ar" || rn == "Cr" || rn == "Gr" || rn == "Ur" || rn == "Tr" ||
       rn == "Ad" || rn == "Cd" || rn == "Gd" || rn == "Td" || rn == "DI")
      return 1;
   return 0;
}

}} // namespace coot::util