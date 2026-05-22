#pragma once
#include <gemmi/model.hpp>
#include "coot/coot-utils/coot-coord-utils.hh"

namespace coot {
namespace util {

// Get the second residue in a chain (index 1).
// Returns nullptr if chain is null or has fewer than 2 residues.
inline const gemmi::Residue* get_second_residue_in_chain_gemmi(const gemmi::Chain& chain) {
    if (chain.residues.size() > 1) {
        return &chain.residues[1];
    }
    return nullptr;
}

inline const gemmi::Residue* get_second_residue_in_chain_gemmi(const gemmi::Model* model, const std::string& chain_id) {
    if (!model) return nullptr;
    for (const auto& chain : model->chains) {
        if (chain.name == chain_id) {
            return get_second_residue_in_chain_gemmi(chain);
        }
    }
    return nullptr;
}

} // namespace util
} // namespace coot