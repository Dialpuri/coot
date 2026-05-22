#pragma once
#include <gemmi/model.hpp>

namespace coot { namespace util {

inline gemmi::CRA next_residue_gemmi(gemmi::CRA cra) {
    if (!cra.residue || !cra.chain)
        return gemmi::CRA{nullptr, nullptr, nullptr};

    auto& chain = *cra.chain;
    auto& residues = chain.residues;
    
    for (size_t i = 0; i < residues.size(); ++i) {
        if (&residues[i] == cra.residue) {
            if (i + 1 < residues.size()) {
                return gemmi::CRA{&chain, &residues[i + 1], nullptr};
            }
            break;
        }
    }
    
    return gemmi::CRA{nullptr, nullptr, nullptr};
}

}} // namespace coot::util