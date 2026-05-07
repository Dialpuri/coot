#pragma once

#include <string>
#include <vector>
#include <set>
#include <gemmi/model.hpp>

namespace coot {

// Collect all unique residue names across all models and chains,
// returned in sorted order (mirroring the original MMDB std::set behavior).
inline std::vector<std::string>
get_types_in_molecule_gemmi(const gemmi::Structure& st) {
    std::set<std::string> s;
    for (const auto& model : st.models) {
        for (const auto& chain : model.chains) {
            for (const auto& residue : chain.residues) {
                s.insert(residue.name);
            }
        }
    }
    return std::vector<std::string>(s.begin(), s.end());
}

} // namespace coot