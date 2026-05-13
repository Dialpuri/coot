#pragma once

#include <utility>
#include <gemmi/model.hpp>

namespace coot { namespace util {

// Return (min_seqnum, max_seqnum) for all residues in the chain.
// For an empty or null chain, returns (9999, -9999) to match MMDB behavior.
inline std::pair<int, int>
min_and_max_residues_gemmi(const gemmi::Chain& chain) {
  std::pair<int, int> p(9999, -9999);
  for (const gemmi::Residue& res : chain.residues) {
    int rn = res.seqid.num.value;
    if (rn < p.first)  p.first  = rn;
    if (rn > p.second) p.second = rn;
  }
  return p;
}

// Overload for nullptr (MMDB compatibility)
inline std::pair<int, int>
min_and_max_residues_gemmi(const gemmi::Chain* chain) {
  if (chain)
    return min_and_max_residues_gemmi(*chain);
  return std::make_pair(9999, -9999);
}

}} // namespace coot::util