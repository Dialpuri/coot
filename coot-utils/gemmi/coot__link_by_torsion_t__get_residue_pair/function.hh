#pragma once

#include <utility>
#include <gemmi/model.hpp>

namespace coot {

class link_by_torsion_t {
public:
    // Original MMDB version:
    // std::pair<mmdb::Residue *, mmdb::Residue *>
    // get_residue_pair(mmdb::Manager *mol);
    
    // Gemmi version - takes a Structure and returns pair of Residue pointers
    std::pair<gemmi::Residue *, gemmi::Residue *>
    get_residue_pair_gemmi(gemmi::Structure& st) {
        std::pair<gemmi::Residue *, gemmi::Residue *> r{nullptr, nullptr};
        
        for (gemmi::Model& model : st.models) {
            for (gemmi::Chain& chain : model.chains) {
                for (gemmi::Residue& residue : chain.residues) {
                    if (r.first) {
                        r.second = &residue;
                        break;
                    } else {
                        r.first = &residue;
                    }
                }
                if (r.first && r.second)
                    break;
            }
            if (r.first && r.second)
                break;
        }
        return r;
    }
};

} // namespace coot