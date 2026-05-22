#pragma once

#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>

namespace coot {
namespace util {

// Return the first residue in the structure, or nullptr if not found
// Returns a CRA (Chain/Residue/Atom) to provide parent context
inline gemmi::CRA get_first_residue_gemmi(const gemmi::Structure& st) {
    if (st.models.empty()) {
        return gemmi::CRA{nullptr, nullptr, nullptr};
    }
    const gemmi::Model& model = st.models[0];
    for (const gemmi::Chain& chain : model.chains) {
        if (!chain.residues.empty()) {
            return {const_cast<gemmi::Chain*>(&chain), const_cast<gemmi::Residue*>(&chain.residues[0]), nullptr};
        }
    }
    return gemmi::CRA{nullptr, nullptr, nullptr};
}

inline gemmi::CRA get_first_residue_gemmi(const gemmi::Model& model) {
    for (const gemmi::Chain& chain : model.chains) {
        if (!chain.residues.empty()) {
            return {const_cast<gemmi::Chain*>(&chain), const_cast<gemmi::Residue*>(&chain.residues[0]), nullptr};
        }
    }
    return gemmi::CRA{nullptr, nullptr, nullptr};
}

} // namespace util
} // namespace coot