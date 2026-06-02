#pragma once

#include <gemmi/model.hpp>
#include <memory>
#include <string>
#include "mini-mol/mini-mol.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__minimol__residue__make_residue/gemmi/function.hh"

namespace coot {
namespace minimol {

/// gemmi port of coot::minimol::fragment::make_chain()
/// Converts a mini-mol fragment into a gemmi::Chain, populating residues.
inline std::unique_ptr<gemmi::Chain> make_chain_gemmi(coot::minimol::fragment const& frag) {
   auto chain_p = std::make_unique<gemmi::Chain>();
   chain_p->name = frag.fragment_id;

   for (int ires = frag.min_res_no(); ires <= frag.max_residue_number(); ires++) {
      auto residue_p = make_residue_gemmi(frag[ires]);
      if (residue_p) {
         chain_p->residues.push_back(std::move(*residue_p));
      }
   }

   return chain_p;
}

} // namespace minimol
} // namespace coot