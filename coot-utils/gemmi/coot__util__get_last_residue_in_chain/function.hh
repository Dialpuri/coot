#pragma once
#include <gemmi/model.hpp>

namespace coot {
namespace util {

// Return the last residue in a chain, or nullptr if chain is empty/invalid.
inline const gemmi::Residue* get_last_residue_in_chain_gemmi(const gemmi::Chain& chain) {
    if (chain.residues.empty()) {
        return nullptr;
    }
    return &chain.residues.back();
}

} // namespace util
} // namespace coot