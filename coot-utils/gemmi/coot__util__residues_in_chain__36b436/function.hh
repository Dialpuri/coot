#pragma once

#include <vector>
#include <string>
#include <gemmi/model.hpp>

namespace coot { namespace util {

// Return all residues in a chain with the given chain ID.
// Returns empty vector if st is nullptr or chain not found.
inline std::vector<gemmi::Residue*>
residues_in_chain_gemmi(gemmi::Structure* st, const std::string& chain_id) {
    std::vector<gemmi::Residue*> v;
    if (!st || st->models.empty())
        return v;

    for (gemmi::Chain& chain : st->models[0].chains) {
        if (chain.name == chain_id) {
            for (gemmi::Residue& residue : chain.residues) {
                v.push_back(&residue);
            }
        }
    }
    return v;
}

}} // namespace coot::util