#pragma once
#include <gemmi/pdb.hpp>
#include <vector>

namespace coot { namespace util {

// Returns the maximum number of residues in any chain in the structure.
// Returns -1 if structure is empty (no models).
inline int max_number_of_residues_in_chain_gemmi(const gemmi::Structure& st) {
   int max_number_of_residues_in_chain = -1;
   if (!st.models.empty()) { 
      for (const gemmi::Model& model : st.models) {
         if (model.chains.empty()) {
            // Skip empty models
            continue;
         }
         for (const gemmi::Chain& chain : model.chains) {
            int nres = chain.residues.size();
            if (nres > max_number_of_residues_in_chain) {
               max_number_of_residues_in_chain = nres;
            }
         }
      }
   }
   return max_number_of_residues_in_chain;
}

}} // namespace coot::util