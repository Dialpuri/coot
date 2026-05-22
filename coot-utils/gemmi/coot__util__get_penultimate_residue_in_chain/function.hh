#pragma once
#include <gemmi/model.hpp>

namespace coot { namespace util {

// Return the penultimate residue in a chain (nullptr if chain has ≤1 residue)
inline const gemmi::Residue*
get_penultimate_residue_in_chain_gemmi(const gemmi::Chain& chain) {
    if (chain.residues.size() > 1) {
        return &chain.residues[chain.residues.size() - 2];
    }
    return nullptr;
}

}} // namespace coot::util