#pragma once

#include <string>
#include <vector>
#include <set>
#include <gemmi/model.hpp>

namespace coot {

// Collect unique residue names from all models/chains in the structure,
// returned in sorted order (via std::set).
inline std::vector<std::string>
get_types_in_molecule_gemmi(const gemmi::Structure &st) {
    std::vector<std::string> v;
    std::set<std::string> s;
    for (const auto &model : st.models) {
        for (const auto &chain : model.chains) {
            for (const auto &residue : chain.residues) {
                s.insert(residue.name);
            }
        }
    }
    for (const auto &item : s) {
        v.push_back(item);
    }
    return v;
}

} // namespace coot