#include "function.hh"
#include <gemmi/neighbor.hpp>
#include <gemmi/contact.hpp>

namespace coot {

std::map<const gemmi::Residue*, std::set<const gemmi::Residue*>>
residues_near_residues_gemmi(
    const std::vector<std::pair<bool, gemmi::CRA>>& residues_vec,
    gemmi::Structure& st,
    float dist_crit) {

    std::map<const gemmi::Residue*, std::set<const gemmi::Residue*>> result;

    if (dist_crit <= 0.0f) return result;
    if (st.models.empty()) return result;

    gemmi::Model& model = st.models[0];

    // Set up NeighborSearch over the entire model
    gemmi::NeighborSearch ns(model, st.cell, static_cast<double>(dist_crit));
    ns.populate(false); // exclude hydrogens

    // For each input residue, find neighboring residues
    for (const auto& [in_input, cra] : residues_vec) {
        if (!cra.residue) continue;
        if (!in_input) continue;

        const gemmi::Residue* key = cra.residue;

        // For each atom in the key residue, find all nearby atoms
        for (const gemmi::Atom& atom : key->atoms) {
            auto neighbors = ns.find_atoms(atom.pos, '\0', 0.0, static_cast<double>(dist_crit));
            for (auto* mark : neighbors) {
                // Skip symmetry mates (image_idx != 0)
                if (mark->image_idx != 0) continue;

                gemmi::const_CRA ncra = mark->to_cra(model);
                const gemmi::Residue* neighbor_res = ncra.residue;
                if (neighbor_res && neighbor_res != key) {
                    result[key].insert(neighbor_res);
                }
            }
        }

        // Only keep the key if it has actual neighbors (matching MMDB behavior)
        if (result.count(key) && result[key].empty()) {
            result.erase(key);
        }
    }

    return result;
}

} // namespace coot