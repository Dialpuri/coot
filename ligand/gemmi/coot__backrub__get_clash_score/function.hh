#pragma once

#include <gemmi/model.hpp>
#include <coot/mini-mol/mini-mol.hh>
#include <clipper/core/coords.h>
#include <vector>
#include <string>
#include <utility>
#include <cmath>
#include <iostream>

namespace coot {
namespace backrub {

/// gemmi port of coot::backrub::get_clash_score
///
/// Computes a clash score between rotamer atoms (from a minimol molecule)
/// and a set of selected sphere atoms (as CRAs).
///
/// Parameters:
///   a_rotamer           - minimol molecule containing the rotamer atoms
///   atom_selection      - CRAs of atoms to check for clashes against
///   chain_id            - chain ID of the current residue (for skipping)
///   resno_1 / resno_2 / resno_3 - sequence numbers of prev/this/next residues
///   water_interaction_mode - 0 = waters excluded from score (collected instead),
///                            1 = waters included in score
///
/// Returns {clash_score, clashing_water_CRAs}
inline std::pair<float, std::vector<gemmi::CRA>> get_clash_score_gemmi(
    const coot::minimol::molecule &a_rotamer,
    const std::vector<gemmi::CRA> &atom_selection,
    const std::string &chain_id,
    int resno_1, int resno_2, int resno_3,
    int water_interaction_mode) {

    std::vector<gemmi::CRA> clashing_waters;
    float clash_score = 0.0f;
    double dist_crit = 3.2;
    double dist_crit_H_bonder = 2.5;
    double dist_crit_sq = dist_crit * dist_crit;
    double dist_crit_H_bonder_sq = dist_crit_H_bonder * dist_crit_H_bonder;

    // Check if rotamer has any atoms (replaces a_rotamer.get_pos().first < 0.0)
    bool has_atoms = false;
    for (unsigned int ifrag = 0; ifrag < a_rotamer.fragments.size(); ifrag++) {
        for (int ires = a_rotamer[ifrag].min_res_no(); ires <= a_rotamer[ifrag].max_residue_number(); ires++) {
            if (a_rotamer[ifrag][ires].n_atoms() > 0) {
                has_atoms = true;
                break;
            }
        }
        if (has_atoms) break;
    }

    if (!has_atoms) {
        std::cout << "ERROR: clash score: there are no atoms in the residue" << std::endl;
    } else {
        for (size_t i = 0; i < atom_selection.size(); i++) {
            gemmi::CRA at = atom_selection[i];
            if (!at.atom) continue;

            clipper::Coord_orth atom_sel_atom_pos(at.atom->pos.x, at.atom->pos.y, at.atom->pos.z);
            int atom_sel_resno = at.residue ? at.residue->seqid.num.value : 0;
            std::string atom_sel_atom_chain = at.chain ? at.chain->name : "";

            // Pad gemmi element name to match MMDB " C" format
            std::string atom_sel_ele_raw = at.atom->element.name();
            std::string atom_sel_ele;
            if (atom_sel_ele_raw.size() == 1)
                atom_sel_ele = " " + atom_sel_ele_raw;
            else
                atom_sel_ele = atom_sel_ele_raw;

            bool count_it = true;
            if (chain_id == atom_sel_atom_chain) {
                if (atom_sel_resno == resno_1 || atom_sel_resno == resno_2 || atom_sel_resno == resno_3) {
                    count_it = false;
                }
            }

            if (count_it) {
                for (unsigned int ifrag = 0; ifrag < a_rotamer.fragments.size(); ifrag++) {
                    for (int ires = a_rotamer[ifrag].min_res_no(); ires <= a_rotamer[ifrag].max_residue_number(); ires++) {
                        for (unsigned int iat = 0; iat < a_rotamer[ifrag][ires].n_atoms(); iat++) {
                            double dlsq = (a_rotamer[ifrag][ires][iat].pos - atom_sel_atom_pos).lengthsq();
                            if (dlsq <= 0.001)
                                dlsq = 0.001;

                            if (dlsq < dist_crit_sq) {
                                bool is_water = false;
                                if (at.residue && at.residue->name == "HOH")
                                    is_water = true;

                                if (!is_water || water_interaction_mode == 1) {
                                    // We take off 1.0/dist_crit_sq so that the clash score goes to 0 at dist_crit.
                                    float extra_clash_score = 100.0f * (static_cast<float>(1.0 / dlsq) - static_cast<float>(1.0 / dist_crit_sq));

                                    // Check if neither atom is carbon (H-bonder logic)
                                    std::string rotamer_elem = a_rotamer[ifrag][ires][iat].element;
                                    if (atom_sel_ele != " C" && rotamer_elem != " C") {
                                        if (dlsq > dist_crit_H_bonder_sq)
                                            extra_clash_score = 0.0f;
                                        else
                                            extra_clash_score = 100.0f * (static_cast<float>(1.0 / dlsq) - static_cast<float>(1.0 / dist_crit_H_bonder_sq));
                                    }
                                    clash_score += extra_clash_score;
                                } else {
                                    double d = std::sqrt(dlsq);
                                    if (d < dist_crit_H_bonder) {
                                        clashing_waters.push_back(at);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return {clash_score, clashing_waters};
}

} // namespace backrub
} // namespace coot