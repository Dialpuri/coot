#pragma once
#include <utility>
#include <gemmi/model.hpp>

namespace coot { namespace util {

// Return (found_residues, min_resno)
inline std::pair<bool, int>
min_resno_in_chain_gemmi(const gemmi::Chain& chain) {
    bool found_residues = false;
    int min_resno = 99999999;
    
    for (const gemmi::Residue& residue : chain.residues) {
        int resno = residue.seqid.num.value;
        if (resno < min_resno) {
            min_resno = resno;
            found_residues = true;
        }
    }
    return std::make_pair(found_residues, min_resno);
}

inline std::pair<bool, int>
min_resno_in_chain_gemmi(const gemmi::Model* model) {
    if (!model || model->chains.empty()) {
        return std::make_pair(false, 99999999);
    }
    int min_resno = 99999999;
    bool found_residues = false;
    for (const gemmi::Chain& chain : model->chains) {
        for (const gemmi::Residue& residue : chain.residues) {
            int resno = residue.seqid.num.value;
            if (resno < min_resno) {
                min_resno = resno;
                found_residues = true;
            }
        }
    }
    return std::make_pair(found_residues, min_resno);
}

} } // namespace coot::util