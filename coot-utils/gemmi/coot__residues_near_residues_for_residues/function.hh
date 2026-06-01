#pragma once

#include <gemmi/model.hpp>
#include <map>
#include <set>

namespace coot {

std::map<gemmi::Residue*, std::set<gemmi::Residue*>>
residues_near_residues_for_residues_gemmi(
    const std::map<gemmi::Residue*, std::set<gemmi::Residue*>>& all_molecule_map,
    const std::set<gemmi::Residue*>& limit_to_these_residues_vec) {

    std::map<gemmi::Residue*, std::set<gemmi::Residue*>> m;

    for (const auto& kv : all_molecule_map) {
        gemmi::Residue* r = kv.first;
        if (limit_to_these_residues_vec.find(r) != limit_to_these_residues_vec.end()) {
            m[r] = kv.second;
        }
    }
    return m;
}

} // namespace coot