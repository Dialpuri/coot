#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <utility>
#include <cmath>
#include <iostream>

namespace coot {

/// gemmi port of coot::get_clash_score
/// Computes a clash score between rotamer atoms and a set of selected atoms.
/// Returns {score, clashing_water_CRAs}.
///
/// - rotamer_atoms: CRAs pointing to the atoms of the rotamer molecule
/// - atom_selection: CRAs pointing to the atoms to check for clashes against
/// - water_interaction_mode: 0 = waters excluded from score (collected instead),
///                           1 = waters included in score
inline std::pair<float, std::vector<gemmi::CRA>> get_clash_score_gemmi(
    const std::vector<gemmi::CRA>& rotamer_atoms,
    const std::vector<gemmi::CRA>& atom_selection,
    int water_interaction_mode) {

    float score = 0.0f;
    float dist_crit = 2.1f;
    std::vector<gemmi::CRA> clashing_waters;

    // Compute mean position and max deviation of rotamer atoms
    // (replaces a_rotamer.get_pos() from minimol)
    gemmi::Vec3 mean_residue_pos = {0, 0, 0};
    size_t n_atoms = 0;
    for (const auto& cra : rotamer_atoms) {
        if (cra.atom) {
            mean_residue_pos += gemmi::Vec3(cra.atom->pos);
            n_atoms++;
        }
    }
    if (n_atoms > 0) {
        mean_residue_pos /= static_cast<double>(n_atoms);
    }

    double max_dev_residue_pos = 0.0;
    for (const auto& cra : rotamer_atoms) {
        if (cra.atom) {
            double d = (gemmi::Vec3(cra.atom->pos) - mean_residue_pos).length();
            if (d > max_dev_residue_pos)
                max_dev_residue_pos = d;
        }
    }

    if (max_dev_residue_pos < 0.0) {
        std::cout << "ERROR: clash score: there are no atoms in the residue" << std::endl;
    } else {
        for (const auto& sel_cra : atom_selection) {
            if (!sel_cra.atom) continue;
            if (!sel_cra.chain) continue;
            if (!sel_cra.residue) continue;

            gemmi::Vec3 sel_atom_pos = sel_cra.atom->pos;
            double d = (sel_atom_pos - mean_residue_pos).length();

            if (d < (max_dev_residue_pos + dist_crit)) {
                int atom_sel_atom_resno = sel_cra.residue->seqid.num.value;
                std::string atom_sel_atom_chain = sel_cra.chain->name;

                for (const auto& rot_cra : rotamer_atoms) {
                    if (!rot_cra.atom) continue;
                    if (!rot_cra.chain) continue;
                    if (!rot_cra.residue) continue;

                    double d_atom = (gemmi::Vec3(rot_cra.atom->pos) - sel_atom_pos).length();
                    std::cout << "  d_atom " << d_atom << "\n";

                    if (d_atom < dist_crit) {
                        int rotamer_resno = rot_cra.residue->seqid.num.value;
                        std::string rotamer_chain = rot_cra.chain->name;

                        std::cout << "comparing rotamer chain :" << rotamer_chain
                                  << ": this res chain " << atom_sel_atom_chain
                                  << " and resnos " << rotamer_resno
                                  << " with this resno " << atom_sel_atom_resno
                                  << std::endl;

                        // Exclude same chain + same residue
                        if (!((rotamer_resno == atom_sel_atom_resno) &&
                              (rotamer_chain == atom_sel_atom_chain))) {
                            // Exclude backbone and hydrogen atoms
                            if (rot_cra.atom->name != " N  " &&
                                rot_cra.atom->name != " C  " &&
                                rot_cra.atom->name != " CA " &&
                                rot_cra.atom->name != " O  " &&
                                rot_cra.atom->name != " H  ") {

                                bool is_water = (sel_cra.residue->name == "HOH");
                                if (!is_water || water_interaction_mode == 1) {
                                    float badness = 100.0f * (1.0f / static_cast<float>(d_atom) - 1.0f / dist_crit);
                                    std::cout << "adding badness " << badness
                                              << " for " << static_cast<const void*>(sel_cra.atom)
                                              << "\n";
                                    if (badness > 100.0f)
                                        badness = 100.0f;
                                    score += badness;
                                } else {
                                    // Collect clashing waters for deletion
                                    clashing_waters.push_back(sel_cra);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return {score, clashing_waters};
}

} // namespace coot