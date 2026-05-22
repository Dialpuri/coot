#pragma once
#include <vector>
#include <string>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>

namespace coot { namespace util {

inline std::vector<gemmi::CRA>
residues_in_molecule_of_type_gemmi(gemmi::Structure& st, const std::string& residue_type) {
   std::vector<gemmi::CRA> v;
   for (gemmi::Model& model : st.models) {
      for (gemmi::Chain& chain : model.chains) {
         for (gemmi::Residue& residue : chain.residues) {
            if (residue.name == residue_type) {
               v.push_back({&chain, &residue, nullptr});
            }
         }
      }
   }
   return v;
}

}} // namespace coot::util