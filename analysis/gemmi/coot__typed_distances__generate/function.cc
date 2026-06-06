#include "function.hh"

#include <gemmi/neighbor.hpp>
#include <gemmi/unitcell.hpp>
#include <cmath>
#include <iostream>

namespace coot {
namespace typed_distances {

void generate_gemmi(
    gemmi::Model& model,
    const std::vector<gemmi::CRA>& selected_atoms,
    std::map<gemmi::CRA, std::map<int, std::vector<float>>, cra_less>& residue_distances_map,
    std::map<gemmi::CRA, std::vector<gemmi::CRA>, cra_less>& residues_within_window
) {
    float max_dist = 5.7;

    if (selected_atoms.empty()) {
        return;
    }

    // Build neighbor search - need non-const Model& and UnitCell
    gemmi::UnitCell cell; // default (non-periodic)
    gemmi::NeighborSearch ns(model, cell, max_dist);
    ns.populate(false); // include_h = false

    // Find residues within window
    find_residues_within_window_gemmi(model, 4, residues_within_window);

    for (int i = 0; i < static_cast<int>(selected_atoms.size()); i++) {
        gemmi::CRA cra1 = selected_atoms[i];
        if (!cra1.atom || !cra1.residue || !cra1.chain) continue;

        auto hits = ns.find_atoms(cra1.atom->pos, '\0', 0.01, max_dist);

        for (auto* mark : hits) {
            gemmi::CRA cra2 = mark->to_cra(model);
            if (!cra2.atom || !cra2.residue || !cra2.chain) continue;

            gemmi::Atom* central_at = cra1.atom;
            gemmi::Atom* neighb_at = cra2.atom;

            // in_self_or_bonded_residue check
            // Same residue
            if (cra1.residue == cra2.residue) continue;

            // Bonded residues: consecutive residues in the same chain
            bool bonded = false;
            if (cra1.chain == cra2.chain) {
                int seq1 = cra1.residue->seqid.num.value;
                int seq2 = cra2.residue->seqid.num.value;
                if (std::abs(seq1 - seq2) <= 1) {
                    bonded = true;
                }
            }
            if (bonded) continue;

            atom_type_t tc = get_type_gemmi(central_at);
            atom_type_t tn = get_type_gemmi(neighb_at);

            double dist_d = distance_gemmi(central_at, neighb_at);
            float dist = static_cast<float>(dist_d);

            int atom_pair_bin_id = get_atom_pair_bin_id(tc, tn);

            if (atom_pair_bin_id == -1) {
                // bad atom-pair-bin-id - skip silently
            } else {
                // Use cra1 as the key for the residue (central residue)
                gemmi::CRA res_cra{cra1.chain, cra1.residue, nullptr};

                residue_distances_map[res_cra][atom_pair_bin_id].push_back(dist);
            }
        }
    }
}

}  // namespace typed_distances
}  // namespace coot