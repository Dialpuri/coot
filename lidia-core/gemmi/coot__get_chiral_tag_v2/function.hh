#pragma once

#include <gemmi/model.hpp>
#include <coot/geometry/protein-geometry.hh>
#include <coot/lidia-core/use-rdkit.hh>
#include <algorithm>
#include <vector>
#include <iostream>
#include <string>

namespace coot {

// Trim whitespace from a 4-char padded atom name to match gemmi's unpadded form
static inline std::string trim_atom_name(const std::string &s) {
    std::string out;
    for (char c : s) {
        if (c != ' ') out += c;
    }
    return out;
}

// Local version of chiral_neighbour_info_t using gemmi::Atom* instead of mmdb::Atom*
struct chiral_neighbour_info_t_gemmi {
    const gemmi::Atom *at;
    int idx_mmcif;
    int idx_atom_list;

    chiral_neighbour_info_t_gemmi(const gemmi::Atom *at_in, int idx_mmcif_in, int idx_atom_list_in)
        : at(at_in), idx_mmcif(idx_mmcif_in), idx_atom_list(idx_atom_list_in) {}

    static bool neighbour_sorter(const chiral_neighbour_info_t_gemmi &a, const chiral_neighbour_info_t_gemmi &b) {
        return a.idx_mmcif < b.idx_mmcif;
    }
};

RDKit::Atom::ChiralType get_chiral_tag_v2_gemmi(
    gemmi::Residue *residue_p,
    const dictionary_residue_restraints_t &restraints,
    const gemmi::Atom *atom_p) {

    RDKit::Atom::ChiralType chiral_tag = RDKit::Atom::CHI_UNSPECIFIED;

    std::string atom_name = atom_p->name;

    std::cout << "Called get_chiral_tag_v2_gemmi() with atom name " << atom_name << std::endl;

    for (unsigned int ich = 0; ich < restraints.chiral_restraint.size(); ich++) {
        const dict_chiral_restraint_t &cr = restraints.chiral_restraint[ich];

        if (trim_atom_name(cr.atom_id_c_4c()) == atom_name) {
            std::vector<chiral_neighbour_info_t_gemmi> neighbs;

            for (int iat = 0; iat < static_cast<int>(residue_p->atoms.size()); iat++) {
                const gemmi::Atom &at = residue_p->atoms[iat];
                std::string atom_name_local(at.name);

                if (trim_atom_name(cr.atom_id_1_4c()) == atom_name_local)
                    neighbs.push_back(chiral_neighbour_info_t_gemmi(&at, iat, 1));
                if (trim_atom_name(cr.atom_id_2_4c()) == atom_name_local)
                    neighbs.push_back(chiral_neighbour_info_t_gemmi(&at, iat, 2));
                if (trim_atom_name(cr.atom_id_3_4c()) == atom_name_local)
                    neighbs.push_back(chiral_neighbour_info_t_gemmi(&at, iat, 3));
            }

            if (neighbs.size() != 3) {
                std::cout << "Errg. Not all chiral neighbours found " << neighbs.size() << " "
                          << atom_name << std::endl;
            } else {
                for (unsigned int ib = 0; ib < restraints.bond_restraint.size(); ib++) {
                    std::string other_atom;
                    const dict_bond_restraint_t &br = restraints.bond_restraint[ib];
                    if (trim_atom_name(br.atom_id_1_4c()) == atom_name)
                        other_atom = trim_atom_name(br.atom_id_2_4c());
                    if (trim_atom_name(br.atom_id_2_4c()) == atom_name)
                        other_atom = trim_atom_name(br.atom_id_1_4c());

                    if (!other_atom.empty()) {
                        for (int iat = 0; iat < static_cast<int>(residue_p->atoms.size()); iat++) {
                            const gemmi::Atom &at = residue_p->atoms[iat];
                            std::string atom_name_local(at.name);

                            if (atom_name_local == other_atom) {
                                bool found = false;
                                for (const auto &n : neighbs) {
                                    if (n.at == &at) {
                                        found = true;
                                        break;
                                    }
                                }

                                if (!found) {
                                    std::cout << atom_name_local << " was not found in neighbs vec" << std::endl;
                                    neighbs.push_back(chiral_neighbour_info_t_gemmi(&at, iat, 0));
                                    std::cout << "neighbs now of size() " << neighbs.size() << std::endl;
                                    break;
                                } else {
                                    std::cout << atom_name_local << " was already in in neighbs vec" << std::endl;
                                }
                            }
                        }
                    }
                }

                if (neighbs.size() != 4) {
                    std::cout << "WARNING:: Errgh. Not we don't have 4 chiral-centre neighbours "
                              << neighbs.size() << std::endl;
                } else {
                    std::sort(neighbs.begin(), neighbs.end(),
                              chiral_neighbour_info_t_gemmi::neighbour_sorter);
                    std::vector<chiral_neighbour_info_t_gemmi> back_neighbs;
                    back_neighbs.push_back(neighbs[1]);
                    back_neighbs.push_back(neighbs[2]);
                    back_neighbs.push_back(neighbs[3]);

                    std::cout << "back_neighbs: (sorted) "
                              << back_neighbs[0].idx_mmcif << " "
                              << back_neighbs[1].idx_mmcif << " "
                              << back_neighbs[2].idx_mmcif << " "
                              << std::endl;
                    std::cout << "back_neighbs:          "
                              << back_neighbs[0].idx_atom_list << " "
                              << back_neighbs[1].idx_atom_list << " "
                              << back_neighbs[2].idx_atom_list << " "
                              << std::endl;

                    bool atom_orders_match = false;
                    // 2 1 0
                    if ((back_neighbs[2].idx_atom_list > back_neighbs[1].idx_atom_list) &&
                        (back_neighbs[1].idx_atom_list > back_neighbs[0].idx_atom_list)) {
                        atom_orders_match = true;
                    }
                    // 0 2 1
                    if ((back_neighbs[0].idx_atom_list > back_neighbs[2].idx_atom_list) &&
                        (back_neighbs[2].idx_atom_list > back_neighbs[1].idx_atom_list)) {
                        atom_orders_match = true;
                    }
                    // 1 0 2
                    if ((back_neighbs[1].idx_atom_list > back_neighbs[0].idx_atom_list) &&
                        (back_neighbs[0].idx_atom_list > back_neighbs[2].idx_atom_list)) {
                        atom_orders_match = true;
                    }

                    if (atom_orders_match) {
                        if (cr.volume_sign == 1)
                            chiral_tag = RDKit::Atom::CHI_TETRAHEDRAL_CW;
                        else
                            chiral_tag = RDKit::Atom::CHI_TETRAHEDRAL_CCW;
                    } else {
                        if (cr.volume_sign == -1)
                            chiral_tag = RDKit::Atom::CHI_TETRAHEDRAL_CW;
                        else
                            chiral_tag = RDKit::Atom::CHI_TETRAHEDRAL_CCW;
                    }
                }
            }
        }
    }

    return chiral_tag;
}

} // namespace coot