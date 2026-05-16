#pragma once

#include <algorithm>
#include <gemmi/model.hpp>

namespace coot {

/// Sort chains in each model of a gemmi::Structure by chain name (ascending).
/// Equivalent to MMDB's SortChains(SORT_CHAIN_ChainID_Asc).
/// PDBCleanup and FinishStructEdit have no gemmi equivalent — gemmi uses
/// immutable structures, so serial/indices are handled on output.
inline void sort_chains_gemmi(gemmi::Structure &st) {
    for (auto &model : st.models) {
        std::sort(model.chains.begin(), model.chains.end(),
                  [](const gemmi::Chain &a, const gemmi::Chain &b) {
                      return a.name < b.name;
                  });
    }
}

} // namespace coot