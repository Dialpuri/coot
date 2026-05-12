#pragma once
#include <vector>
#include <string>
#include <gemmi/model.hpp>

namespace coot {

// Merge chains from other structures into the first structure's first model
// Only copies the first chain from each structure in st_others
void merge_molecules_gemmi(gemmi::Structure& st_first, std::vector<gemmi::Structure> st_others) {
    // Helper to convert chain index to chain ID (A0, A1, ..., A9, B0, ...)
    auto index_to_chain_id = [](unsigned int idx) {
        unsigned int tens = idx / 10;
        unsigned int units = idx - 10 * tens;
        std::string abc("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
        char t = abc[tens];
        std::string result = std::string(1, t) + std::to_string(units);
        return result;
    };

    if (st_first.models.empty())
        return;

    gemmi::Model& model_first = st_first.models[0];
    unsigned int chain_idx = 0;

    for (const auto& st_other : st_others) {
        if (st_other.models.empty())
            continue;

        const gemmi::Model& model_other = st_other.models[0];
        for (const gemmi::Chain& chain : model_other.chains) {
            gemmi::Chain copy_chain = chain;
            std::string ch_id = index_to_chain_id(chain_idx);
            copy_chain.name = ch_id;
            model_first.chains.push_back(copy_chain);
            chain_idx++;
        }
    }
}

} // namespace coot