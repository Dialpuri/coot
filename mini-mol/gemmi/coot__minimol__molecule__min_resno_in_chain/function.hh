#pragma once
#include <utility>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>

namespace coot { namespace minimol { namespace molecule {

// Original: mmdb::Chain* → std::pair<bool, int>
// Ported: gemmi::Chain& → std::pair<bool, int>
// Returns (found, min_resno)
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

}}}