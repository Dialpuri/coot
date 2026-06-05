#include "function.hh"

#include <iostream>

namespace coot {

int add_link_angle_gemmi(std::string link_type,
                         gemmi::CRA first, gemmi::CRA second,
                         short int is_fixed_first,
                         short int is_fixed_second,
                         const coot::protein_geometry &geom) {

    int nangle = 0;

    // Collect atoms from first and second residues
    std::vector<gemmi::Atom*> first_sel;
    std::vector<gemmi::Atom*> second_sel;

    if (first.residue) {
        for (auto &atom : first.residue->atoms) {
            first_sel.push_back(&atom);
        }
    }
    if (second.residue) {
        for (auto &atom : second.residue->atoms) {
            second_sel.push_back(&atom);
        }
    }

    int n_first_res_atoms = static_cast<int>(first_sel.size());
    int n_second_res_atoms = static_cast<int>(second_sel.size());

    if (n_first_res_atoms <= 0) {
        std::cout << "no atoms in first residue!? " << std::endl;
    }
    if (n_second_res_atoms <= 0) {
        std::cout << "no atoms in second residue!? " << std::endl;
    }

    for (int i = 0; i < static_cast<int>(geom.link_size()); i++) {
        if (geom.link(i).link_id == link_type) {
            for (unsigned int j = 0; j < geom.link(i).link_angle_restraint.size(); j++) {

                std::string dict_atom_1 = geom.link(i).link_angle_restraint[j].atom_id_1_4c();
                std::string dict_atom_2 = geom.link(i).link_angle_restraint[j].atom_id_2_4c();
                std::string dict_atom_3 = geom.link(i).link_angle_restraint[j].atom_id_3_4c();

                gemmi::Atom** atom_1_sel;
                int n_atom_1;
                gemmi::Atom** atom_2_sel;
                int n_atom_2;
                gemmi::Atom** atom_3_sel;
                int n_atom_3;

                if (geom.link(i).link_angle_restraint[j].atom_1_comp_id == 1) {
                    atom_1_sel = first_sel.data();
                    n_atom_1 = n_first_res_atoms;
                } else {
                    atom_1_sel = second_sel.data();
                    n_atom_1 = n_second_res_atoms;
                }

                if (geom.link(i).link_angle_restraint[j].atom_2_comp_id == 1) {
                    atom_2_sel = first_sel.data();
                    n_atom_2 = n_first_res_atoms;
                } else {
                    atom_2_sel = second_sel.data();
                    n_atom_2 = n_second_res_atoms;
                }

                if (geom.link(i).link_angle_restraint[j].atom_3_comp_id == 1) {
                    atom_3_sel = first_sel.data();
                    n_atom_3 = n_first_res_atoms;
                } else {
                    atom_3_sel = second_sel.data();
                    n_atom_3 = n_second_res_atoms;
                }

                for (int ifat = 0; ifat < n_atom_1; ifat++) {
                    // Remove trailing spaces from gemmi atom name, then expand to 4-char
                    std::string name_1 = atom_1_sel[ifat]->name;
                    while (!name_1.empty() && name_1.back() == ' ') name_1.pop_back();
                    std::string pdb_atom_name_1 = coot::atom_id_mmdb_expand(name_1);
                    
                    if (pdb_atom_name_1 == dict_atom_1) {
                        for (int isat = 0; isat < n_atom_2; isat++) {
                            std::string name_2 = atom_2_sel[isat]->name;
                            while (!name_2.empty() && name_2.back() == ' ') name_2.pop_back();
                            std::string pdb_atom_name_2 = coot::atom_id_mmdb_expand(name_2);

                            if (pdb_atom_name_2 == dict_atom_2) {
                                for (int itat = 0; itat < n_atom_3; itat++) {
                                    std::string name_3 = atom_3_sel[itat]->name;
                                    while (!name_3.empty() && name_3.back() == ' ') name_3.pop_back();
                                    std::string pdb_atom_name_3 = coot::atom_id_mmdb_expand(name_3);

                                    if (pdb_atom_name_3 == dict_atom_3) {

                                        char alt_conf_1 = atom_1_sel[ifat]->altloc;
                                        char alt_conf_2 = atom_2_sel[isat]->altloc;
                                        char alt_conf_3 = atom_3_sel[itat]->altloc;

                                        // Normalization: gemmi uses ' ' for "no altloc", same as MMDB's ""
                                        if (((alt_conf_1 == alt_conf_2) && (alt_conf_1 == alt_conf_3)) ||
                                            ((alt_conf_1 == alt_conf_2) && (alt_conf_3 == ' ')) ||
                                            ((alt_conf_2 == alt_conf_3) && (alt_conf_1 == ' ')) ||
                                            ((alt_conf_1 == ' ') && (alt_conf_2 == ' ')) ||
                                            ((alt_conf_2 == ' ') && (alt_conf_3 == ' '))) {

                                            nangle++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return nangle;
}

} // namespace coot