#pragma once

#include <string>
#include <vector>
#include <gemmi/model.hpp>

namespace coot {
namespace molecule_t {

// Select residues in chain_id within the range [resno_start, resno_end] (inclusive)
// Returns vector of CRA (Chain-Residue-Atom) for each matching residue
inline std::vector<gemmi::CRA>
select_residues_gemmi(const std::string &chain_id, int resno_start, int resno_end, gemmi::Structure &st) {
   std::vector<gemmi::CRA> rv;

   if (st.models.empty()) {
      return rv;
   }

   gemmi::Model &model = st.models[0];
   for (auto &chain : model.chains) {
      if (chain.name == chain_id) {
         for (auto &residue : chain.residues) {
            int res_no_this = residue.seqid.num.value;
            if (res_no_this >= resno_start && res_no_this <= resno_end) {
               rv.push_back({&chain, &residue, nullptr});
            }
         }
      }
   }
   return rv;
}

} // namespace molecule_t
} // namespace coot