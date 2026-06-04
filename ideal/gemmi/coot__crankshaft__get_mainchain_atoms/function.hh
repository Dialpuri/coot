#pragma once

#include <gemmi/model.hpp>
#include <vector>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__is_main_chain_p/gemmi/function.hh"

namespace coot {
namespace crankshaft {

// gemmi port of coot::crankshaft::get_mainchain_atoms
// Returns pointers to main-chain atoms from res_1 (only if both residues have atoms)
inline std::vector<const gemmi::Atom*> get_mainchain_atoms_gemmi(const gemmi::Residue& res_1, const gemmi::Residue& res_2) {
  std::vector<const gemmi::Atom*> v;

  if (!res_1.atoms.empty() && !res_2.atoms.empty()) {
    for (const gemmi::Atom& at : res_1.atoms) {
      if (coot::is_main_chain_p_gemmi(at, res_1)) {
        v.push_back(&at);
      }
    }
  }

  return v;
}

} // namespace crankshaft
} // namespace coot