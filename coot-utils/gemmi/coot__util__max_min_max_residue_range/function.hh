#pragma once

#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>

namespace coot { namespace util {

int max_min_max_residue_range_gemmi(const gemmi::Structure& st) {
    int max_min_max = -1;
    
    for (const gemmi::Model& model : st.models) {
        for (const gemmi::Chain& chain : model.chains) {
            if (chain.residues.empty()) {
                continue;
            }
            
            int min_resno = 99999;
            int max_resno = -99999;
            
            for (const gemmi::Residue& residue : chain.residues) {
                int this_resno = residue.seqid.num.value;
                if (this_resno > max_resno) {
                    max_resno = this_resno;
                }
                if (this_resno < min_resno) {
                    min_resno = this_resno;
                }
            }
            
            int range = max_resno - min_resno + 1;
            if (range > max_min_max) {
                max_min_max = range;
            }
        }
    }
    
    return max_min_max;
}

}} // namespace coot::util