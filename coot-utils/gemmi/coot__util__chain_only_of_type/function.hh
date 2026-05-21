#pragma once
#include <gemmi/model.hpp>
#include <string>

namespace coot { namespace util {

// Return the first chain in the structure where ALL residues have the same name.
// Returns nullptr if no such chain exists or if the structure has no models.
inline const gemmi::Chain* chain_only_of_type_gemmi(const gemmi::Structure& st, const std::string& residue_type) {
    if (st.models.empty())
        return nullptr;

    for (const gemmi::Model& model : st.models) {
        for (const gemmi::Chain& chain : model.chains) {
            bool all_same_type = true;
            for (const gemmi::Residue& res : chain.residues) {
                if (res.name != residue_type) {
                    all_same_type = false;
                    break;
                }
            }
            if (all_same_type && !chain.residues.empty()) {
                return &chain;
            }
        }
    }
    return nullptr;
}

}} // namespace coot::util