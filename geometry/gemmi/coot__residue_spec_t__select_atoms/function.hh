#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <string>

#include "geometry/residue-and-atom-specs.hh"

namespace coot {

inline std::vector<gemmi::Atom*> select_atoms_gemmi(const residue_spec_t& spec, gemmi::Structure* st) {
  std::vector<gemmi::Atom*> atoms;
  if (st == nullptr) return atoms;

  // Normalize insertion code: MMDB uses "" for "no inscode", gemmi uses ' '
  auto norm = [](const std::string& ic) { return ic.empty() ? std::string(" ") : ic; };
  std::string norm_ins = norm(spec.ins_code);

  for (gemmi::Model& model : st->models) {
    for (gemmi::Chain& chain : model.chains) {
      if (chain.name != spec.chain_id) continue;
      for (gemmi::Residue& residue : chain.residues) {
        if (residue.seqid.num.value != spec.res_no) continue;
        if (std::string(1, residue.seqid.icode) != norm_ins) continue;
        for (gemmi::Atom& atom : residue.atoms) {
          atoms.push_back(&atom);
        }
        break; // Found the residue, move to next chain/model
      }
    }
  }
  return atoms;
}

} // namespace coot