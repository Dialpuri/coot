#pragma once
#include <vector>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>

namespace coot { namespace util {

std::vector<gemmi::Residue*>
residues_with_insertion_codes_gemmi(const gemmi::Structure& st) {
    std::vector<gemmi::Residue*> v;
    
    if (st.models.empty()) return v;
    
    gemmi::Model& model = const_cast<gemmi::Model&>(st.models[0]);
    for (gemmi::Chain& chain : model.chains) {
        for (gemmi::Residue& residue : chain.residues) {
            if (residue.seqid.icode != ' ') {
                v.push_back(&residue);
            }
        }
    }
    
    return v;
}

}} // namespace coot::util