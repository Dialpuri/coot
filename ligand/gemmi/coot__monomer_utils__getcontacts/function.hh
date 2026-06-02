#pragma once

#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <gemmi/unitcell.hpp>
#include <vector>
#include <utility>
#include <unordered_map>

namespace coot {

namespace monomer_utils {

// Port of coot::monomer_utils::getcontacts to gemmi.
//
// Finds contacts (atom pairs within 0.1–1.9 Å) among selected atoms,
// using gemmi::NeighborSearch for spatial indexing.
// Returns vector of (id1, id2) pairs where ids are indices into selected_atoms.
std::vector<std::pair<int, int>> getcontacts_gemmi(
    gemmi::Model& model,
    std::vector<gemmi::Atom*>& selected_atoms)
{
    std::vector<std::pair<int, int>> contacts;
    double min_dist = 0.1;
    double max_dist = 1.9;

    size_t n = selected_atoms.size();
    if (n < 2)
        return contacts;

    // Map atom pointer → index in selected_atoms for fast lookup
    std::unordered_map<gemmi::Atom*, int> atom_to_idx;
    for (size_t i = 0; i < n; ++i) {
        atom_to_idx[selected_atoms[i]] = static_cast<int>(i);
    }

    // Use NeighborSearch on the model (default UnitCell = non-crystallographic)
    gemmi::UnitCell cell;
    gemmi::NeighborSearch ns(model, cell, max_dist);
    ns.populate(true); // include H atoms

    for (size_t i = 0; i < n; ++i) {
        auto hits = ns.find_atoms(selected_atoms[i]->pos, '\0', min_dist, max_dist);
        for (auto* mark : hits) {
            gemmi::CRA cra = mark->to_cra(model);
            if (cra.atom != nullptr && cra.atom != selected_atoms[i]) {
                auto it = atom_to_idx.find(cra.atom);
                if (it != atom_to_idx.end()) {
                    contacts.emplace_back(static_cast<int>(i), it->second);
                }
            }
        }
    }

    return contacts;
}

} // namespace monomer_utils
} // namespace coot