#pragma once

#include <vector>
#include <gemmi/model.hpp>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__like_a_helix/gemmi/function.hh"

namespace coot::daca {

void
fill_helix_flags_gemmi(gemmi::Model *model_p, std::vector<gemmi::CRA> &helical_residues) {

   if (!model_p) return;

   for (gemmi::Chain &chain : model_p->chains) {
      // Collect all residues from this chain
      std::vector<gemmi::Residue *> selected_residues;
      for (gemmi::Residue &res : chain.residues) {
         selected_residues.push_back(&res);
      }

      std::vector<gemmi::Residue *> helical_residues_in_chain =
          coot::like_a_helix_gemmi(selected_residues);

      for (gemmi::Residue *res : helical_residues_in_chain) {
         gemmi::CRA cra{&chain, res, nullptr};
         helical_residues.push_back(cra);
      }
   }
}

} // namespace coot::daca