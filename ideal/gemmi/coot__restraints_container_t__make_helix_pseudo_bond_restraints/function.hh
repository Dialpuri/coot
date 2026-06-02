#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <iostream>

namespace coot {

enum restraint_type_enum { BOND_RESTRAINT = 0 };

class restraints_container_t {
public:
    bool from_residue_vector;
    std::string chain_id_save;
    int istart_res;
    int iend_res;

    gemmi::Model* model;

    // residue storage for from_residue_vector mode
    std::vector<gemmi::Residue*> residues_;

    struct simple_restraint {
        restraint_type_enum type;
        int index1;
        int index2;
        std::vector<bool> fixed_flags;
        double ideal_value;
        double esd;
        double weight;
    };
    std::vector<simple_restraint> restraints_;

    restraints_container_t(gemmi::Model* m = nullptr)
        : from_residue_vector(true),
          chain_id_save(""),
          istart_res(0),
          iend_res(0),
          model(m) {}

    int size() const { return static_cast<int>(restraints_.size()); }

    std::vector<bool> make_fixed_flags(int, int) const {
        return std::vector<bool>{};
    }

    void add(restraint_type_enum type,
             int i1,
             int i2,
             std::vector<bool> fixed_flags,
             double ideal,
             double esd,
             double weight) {
        simple_restraint r;
        r.type = type;
        r.index1 = i1;
        r.index2 = i2;
        r.fixed_flags = fixed_flags;
        r.ideal_value = ideal;
        r.esd = esd;
        r.weight = weight;
        restraints_.push_back(r);
    }

    void make_helix_pseudo_bond_restraints_from_res_vec() {
        // When residue vector is empty, no restraints are added.
        // Mirrors the original: iterates stored residues_, skipping alts,
        // and builds N->O(i-4) and N->O(i-3) helix bonds.
        float pseudo_bond_esd = 0.04;

        std::vector<gemmi::Residue*> sel_residues;
        for (gemmi::Residue* res : residues_) {
            // Skip residues with alt_loc
            bool has_alt = false;
            for (const gemmi::Atom& atom : res->atoms) {
                if (atom.altloc != '\0' && atom.altloc != ' ') {
                    has_alt = true;
                    break;
                }
            }
            if (!has_alt) {
                sel_residues.push_back(res);
            }
        }

        if (sel_residues.size() == 0) return;

        for (size_t i = 4; i < sel_residues.size(); i++) {
            gemmi::Residue* res1 = sel_residues[i];
            for (const gemmi::Atom& atom1 : res1->atoms) {
                if (atom1.name != "N") continue;

                int index1 = static_cast<int>(i);
                int index2 = -1;

                // N(i) -> O(i-4) distance 2.91
                gemmi::Residue* contact_res = sel_residues[i - 4];
                if (res1->seqid.num.value == (contact_res->seqid.num.value + 4)) {
                    for (const gemmi::Atom& atom2 : contact_res->atoms) {
                        if (atom2.name == "O") {
                            index2 = static_cast<int>(i - 4);
                            std::vector<bool> fixed_flags =
                                make_fixed_flags(index1, index2);
                            add(BOND_RESTRAINT, index1, index2, fixed_flags,
                                2.91, pseudo_bond_esd, 1.2);
                            std::cout << "Helix Bond restraint (" << atom1.name
                                      << " " << res1->seqid.num.value
                                      << ") to (" << atom2.name << " "
                                      << contact_res->seqid.num.value
                                      << ") 2.91" << std::endl;
                        }
                    }
                }

                // N(i) -> O(i-3) distance 3.18
                contact_res = sel_residues[i - 3];
                if (res1->seqid.num.value == (contact_res->seqid.num.value + 3)) {
                    for (const gemmi::Atom& atom2 : contact_res->atoms) {
                        if (atom2.name == "O") {
                            index2 = static_cast<int>(i - 3);
                            std::vector<bool> fixed_flags =
                                make_fixed_flags(index1, index2);
                            add(BOND_RESTRAINT, index1, index2, fixed_flags,
                                3.18, pseudo_bond_esd, 1.2);
                            std::cout << "Helix Bond restraint (" << atom1.name
                                      << " " << res1->seqid.num.value
                                      << ") to (" << atom2.name << " "
                                      << contact_res->seqid.num.value
                                      << ") 3.18" << std::endl;
                        }
                    }
                }
            }
        }
    }

    void make_helix_pseudo_bond_restraints_gemmi() {
        // somewhat hacky
        if (from_residue_vector) {
            make_helix_pseudo_bond_restraints_from_res_vec();
            return;
        }

        float pseudo_bond_esd = 0.04;

        if (!model) return;

        for (gemmi::Chain& chain : model->chains) {
            if (chain.name != chain_id_save) continue;

            // Collect residues in the sequence range
            // Skip residues with alt_loc (mirrors original: altLocs="")
            std::vector<gemmi::Residue*> sel_residues;
            for (gemmi::Residue& res : chain.residues) {
                int seq_num = res.seqid.num.value;
                if (seq_num >= istart_res && seq_num <= iend_res) {
                    bool has_alt = false;
                    for (const gemmi::Atom& atom : res.atoms) {
                        if (atom.altloc != '\0' && atom.altloc != ' ') {
                            has_alt = true;
                            break;
                        }
                    }
                    if (!has_alt) {
                        sel_residues.push_back(&res);
                    }
                }
            }

            if (sel_residues.size() == 0) continue;

            for (size_t i = 4; i < sel_residues.size(); i++) {
                gemmi::Residue* res1 = sel_residues[i];
                for (const gemmi::Atom& atom1 : res1->atoms) {
                    if (atom1.name != "N") continue;

                    int index1 = static_cast<int>(i);
                    int index2 = -1;

                    // N(i) -> O(i-4) distance 2.91
                    gemmi::Residue* contact_res = sel_residues[i - 4];
                    if (res1->seqid.num.value ==
                        (contact_res->seqid.num.value + 4)) {
                        for (const gemmi::Atom& atom2 : contact_res->atoms) {
                            if (atom2.name == "O") {
                                index2 = static_cast<int>(i - 4);
                                std::vector<bool> fixed_flags =
                                    make_fixed_flags(index1, index2);
                                add(BOND_RESTRAINT, index1, index2, fixed_flags,
                                    2.91, pseudo_bond_esd, 1.2);
                                std::cout << "Helix Bond restraint ("
                                          << atom1.name << " "
                                          << res1->seqid.num.value << ") to ("
                                          << atom2.name << " "
                                          << contact_res->seqid.num.value
                                          << ") 2.91" << std::endl;
                            }
                        }
                    }

                    // N(i) -> O(i-3) distance 3.18
                    contact_res = sel_residues[i - 3];
                    if (res1->seqid.num.value ==
                        (contact_res->seqid.num.value + 3)) {
                        for (const gemmi::Atom& atom2 : contact_res->atoms) {
                            if (atom2.name == "O") {
                                index2 = static_cast<int>(i - 3);
                                std::vector<bool> fixed_flags =
                                    make_fixed_flags(index1, index2);
                                add(BOND_RESTRAINT, index1, index2, fixed_flags,
                                    3.18, pseudo_bond_esd, 1.2);
                                std::cout << "Helix Bond restraint ("
                                          << atom1.name << " "
                                          << res1->seqid.num.value << ") to ("
                                          << atom2.name << " "
                                          << contact_res->seqid.num.value
                                          << ") 3.18" << std::endl;
                            }
                        }
                    }
                }
            }
        }
    }
};

}  // namespace coot