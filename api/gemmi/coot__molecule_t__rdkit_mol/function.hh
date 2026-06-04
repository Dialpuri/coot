#pragma once

#include <string>
#include <gemmi/model.hpp>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__molecule_t__cid_to_residue/gemmi/function.hh"

namespace coot {
namespace molecule_t {

// Port of coot::molecule_t::rdkit_mol — stub that always returns nullptr
// (the original never actually constructs an RDKit molecule)
inline void* rdkit_mol_gemmi(const std::string& ligand_cid, gemmi::Model& model) {
  auto* residue_p = coot::cid_to_residue_gemmi(ligand_cid, model);
  (void)residue_p;  // Stub: original code gets residue but does not use it
  return nullptr;
}

}  // namespace molecule_t
}  // namespace coot