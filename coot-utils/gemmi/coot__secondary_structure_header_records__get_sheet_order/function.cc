#include "function.hh"
#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <set>
#include <algorithm>
#include <cmath>
#include <map>

// Helper: gemmi version of get_strand_sense
// Determines if two strands run parallel or anti-parallel
static coot::secondary_structure_header_records::strand_relation_t::sense_t
get_strand_sense_gemmi(const std::vector<gemmi::Residue*> &strand_1,
                       const std::vector<gemmi::Residue*> &strand_2) {
    coot::secondary_structure_header_records::strand_relation_t::sense_t s =
        coot::secondary_structure_header_records::strand_relation_t::NO_RESULT;

    if (strand_1.size() > 1 && strand_2.size() > 1) {
        const gemmi::Residue *sr10 = strand_1[0];
        const gemmi::Residue *sr1e = strand_1.back();
        const gemmi::Residue *sr20 = strand_2[0];
        const gemmi::Residue *sr2e = strand_2.back();

        // Compute average position of each residue
        auto avg_pos = [](const gemmi::Residue *res) -> gemmi::Vec3 {
            gemmi::Vec3 center(0, 0, 0);
            int count = 0;
            for (const gemmi::Atom &at : res->atoms) {
                center += at.pos;
                count++;
            }
            if (count == 0) return gemmi::Vec3(0, 0, 0);
            return center / count;
        };

        gemmi::Vec3 sr10_pt = avg_pos(sr10);
        gemmi::Vec3 sr1e_pt = avg_pos(sr1e);
        gemmi::Vec3 sr20_pt = avg_pos(sr20);
        gemmi::Vec3 sr2e_pt = avg_pos(sr2e);

        gemmi::Vec3 v1 = sr1e_pt - sr10_pt;
        gemmi::Vec3 v2 = sr2e_pt - sr20_pt;

        double len1 = std::sqrt(v1.length_sq());
        double len2 = std::sqrt(v2.length_sq());

        if (len1 > 1e-6 && len2 > 1e-6) {
            gemmi::Vec3 v1u = v1 / len1;
            gemmi::Vec3 v2u = v2 / len2;

            double cos_theta = v1u.dot(v2u);
            if (cos_theta > 0)
                s = coot::secondary_structure_header_records::strand_relation_t::PARALLEL;
            else
                s = coot::secondary_structure_header_records::strand_relation_t::ANTI_PARALLEL;
        }
    }
    return s;
}

std::vector<std::vector<coot::secondary_structure_header_records::strand_relation_t>>
get_sheet_order_gemmi(
    gemmi::Structure &st,
    gemmi::Model &model,
    const std::vector<std::vector<gemmi::Residue*>> &strands_with_residues) {

    std::vector<std::vector<coot::secondary_structure_header_records::strand_relation_t>> ordered_strands;

    // Collect N and O atoms from all residues in the strands
    std::vector<gemmi::Atom*> N_atoms_vec;
    std::vector<gemmi::Atom*> O_atoms_vec;

    for (size_t i = 0; i < strands_with_residues.size(); i++) {
        for (size_t ii = 0; ii < strands_with_residues[i].size(); ii++) {
            gemmi::Residue *residue_p = strands_with_residues[i][ii];
            for (gemmi::Atom &at : residue_p->atoms) {
                std::string atom_name = at.name;
                // gemmi atom names from PDB are not padded to 4 chars
                if (atom_name == "O" || atom_name == " O" || atom_name == " O " || atom_name == " O  ")
                    O_atoms_vec.push_back(&at);
                if (atom_name == "N" || atom_name == " N" || atom_name == " N " || atom_name == " N  ")
                    N_atoms_vec.push_back(&at);
            }
        }
    }

    if (N_atoms_vec.empty() || O_atoms_vec.empty()) {
        return ordered_strands;
    }

    double max_dist = 3.5;
    double max_dist_sq = max_dist * max_dist;

    // Build strand ordering: for each strand, a set of its neighboring strand relations
    std::vector<std::set<coot::secondary_structure_header_records::strand_relation_t>>
        strand_ordering(strands_with_residues.size());

    // Build a residue pointer to strand index map for quick lookup
    // Also find residue for each N atom by building a map
    std::map<gemmi::Atom*, gemmi::Residue*> atom_to_residue;
    for (size_t i = 0; i < strands_with_residues.size(); i++) {
        for (gemmi::Residue* res : strands_with_residues[i]) {
            for (gemmi::Atom& at : res->atoms) {
                atom_to_residue[&at] = res;
            }
        }
    }

    // For each N atom, find nearby O atoms by distance
    for (gemmi::Atom* n_atom : N_atoms_vec) {
        gemmi::Residue* r_1 = atom_to_residue[n_atom];
        if (!r_1) continue;

        for (gemmi::Atom* o_atom : O_atoms_vec) {
            double dist_sq = (n_atom->pos - o_atom->pos).length_sq();
            if (dist_sq > max_dist_sq) continue;
            if (dist_sq < 0.01) continue; // min_dist 0.1

            gemmi::Residue* r_2 = atom_to_residue[o_atom];
            if (!r_2) continue;
            if (r_1 == r_2) continue;

            int rn_1 = r_1->seqid.num.value;
            int rn_2 = r_2->seqid.num.value;
            if (std::abs(rn_1 - rn_2) <= 2) continue;

            // Find which strands contain these residues
            for (size_t i = 0; i < strands_with_residues.size(); i++) {
                if (std::find(strands_with_residues[i].begin(),
                              strands_with_residues[i].end(), r_1) != strands_with_residues[i].end()) {
                    for (size_t j = 0; j < strands_with_residues.size(); j++) {
                        if (j != i) {
                            if (std::find(strands_with_residues[j].begin(),
                                          strands_with_residues[j].end(), r_2) != strands_with_residues[j].end()) {
                                coot::secondary_structure_header_records::strand_relation_t::sense_t
                                    sense = get_strand_sense_gemmi(strands_with_residues[i], strands_with_residues[j]);
                                coot::secondary_structure_header_records::strand_relation_t sr_1(j, sense);
                                coot::secondary_structure_header_records::strand_relation_t sr_2(i, sense);

                                strand_ordering[i].insert(sr_1);
                                strand_ordering[j].insert(sr_2);
                            }
                        }
                    }
                }
            }
        }
    }

    std::vector<unsigned int> done_strands;

    // Build ordered sheet starting from strands with 1 neighbor
    for (size_t i = 0; i < strand_ordering.size(); i++) {
        std::vector<coot::secondary_structure_header_records::strand_relation_t> sheet_order_local;

        if (strand_ordering[i].size() == 1) { // start of a sheet
            if (std::find(done_strands.begin(), done_strands.end(), i) == done_strands.end()) {
                coot::secondary_structure_header_records::strand_relation_t
                    sr_first(i, coot::secondary_structure_header_records::strand_relation_t::FIRST);
                sheet_order_local.push_back(sr_first);
                done_strands.push_back(i);
                unsigned int this_strand = i;
                bool cont = true;
                while (cont) {
                    const std::set<coot::secondary_structure_header_records::strand_relation_t> &sr =
                        strand_ordering[this_strand];
                    std::set<coot::secondary_structure_header_records::strand_relation_t>::const_iterator it;
                    bool found_one = false;
                    for (it = sr.begin(); it != sr.end(); it++) {
                        if (std::find(sheet_order_local.begin(), sheet_order_local.end(), *it) ==
                            sheet_order_local.end()) {
                            sheet_order_local.push_back(*it);
                            done_strands.push_back(it->strand_idx);
                            this_strand = it->strand_idx;
                            found_one = true;
                            break;
                        }
                    }
                    if (!found_one) {
                        cont = false;
                    }
                }
            }
        }
        if (!sheet_order_local.empty()) {
            ordered_strands.push_back(sheet_order_local);
        }
    }

    return ordered_strands;
}