#pragma once
#include <gemmi/model.hpp>

namespace coot { namespace util {

// Return the first residue in the chain, or nullptr if the chain has no residues.
// Same semantics as the original MMDB function.
inline const gemmi::Residue* get_first_residue_in_chain_gemmi(const gemmi::Chain& chain) {
    if (chain.residues.empty()) {
        return nullptr;
    }
    return &chain.residues.front();
}

}} // namespace coot::util