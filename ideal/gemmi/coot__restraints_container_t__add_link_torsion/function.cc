#include "function.hh"

#include <iostream>
#include <vector>

namespace coot {

int add_link_torsion_gemmi(std::string link_type,
                           gemmi::CRA first,
                           gemmi::CRA second,
                           short int is_fixed_first,
                           short int is_fixed_second,
                           const coot::protein_geometry &geom) {

    int n_torsions = 0;

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

    std::vector<bool> fixed_flag(4);
    fixed_flag[0] = false;
    fixed_flag[1] = false;
    fixed_flag[2] = false;
    fixed_flag[3] = false;

    for (int i = 0; i < static_cast<int>(geom.link_size()); i++) {
        if (geom.link(i).link_id == link_type) {
            for (unsigned int j = 0; j < geom.link(i).link_torsion_restraint.size(); j++) {

                const auto &ltr = geom.link(i).link_torsion_restraint[j];

                // Determine which residue each atom comes from
                gemmi::Atom** atom_1_sel;
                int n_atom_1;
                if (ltr.atom_1_comp_id == 1) {
                    atom_1_sel = first_sel.data();
                    n_atom_1 = n_first_res_atoms;
                    fixed_flag[0] = is_fixed_first;
                } else {
                    atom_1_sel = second_sel.data();
                    n_atom_1 = n_second_res_atoms;
                    fixed_flag[0] = is_fixed_second;
                }

                gemmi::Atom** atom_2_sel;
                int n_atom_2;
                if (ltr.atom_2_comp_id == 1) {
                    atom_2_sel = first_sel.data();
                    n_atom_2 = n_first_res_atoms;
                    fixed_flag[1] = is_fixed_first;
                } else {
                    atom_2_sel = second_sel.data();
                    n_atom_2 = n_second_res_atoms;
                    fixed_flag[1] = is_fixed_second;
                }

                gemmi::Atom** atom_3_sel;
                int n_atom_3;
                if (ltr.atom_3_comp_id == 1) {
                    atom_3_sel = first_sel.data();
                    n_atom_3 = n_first_res_atoms;
                    fixed_flag[2] = is_fixed_first;
                } else {
                    atom_3_sel = second_sel.data();
                    n_atom_3 = n_second_res_atoms;
                    fixed_flag[2] = is_fixed_second;
                }

                gemmi::Atom** atom_4_sel;
                int n_atom_4;
                if (ltr.atom_4_comp_id == 1) {
                    atom_4_sel = first_sel.data();
                    n_atom_4 = n_first_res_atoms;
                    fixed_flag[3] = is_fixed_first;
                } else {
                    atom_4_sel = second_sel.data();
                    n_atom_4 = n_second_res_atoms;
                    fixed_flag[3] = is_fixed_second;
                }

                for (int ifat = 0; ifat < n_atom_1; ifat++) {
                    std::string name_1 = atom_1_sel[ifat]->name;
                    while (!name_1.empty() && name_1.back() == ' ') name_1.pop_back();
                    std::string pdb_atom_name_1 = coot::atom_id_mmdb_expand(name_1);

                    if (pdb_atom_name_1 == ltr.atom_id_1_4c()) {
                        for (int isat = 0; isat < n_atom_2; isat++) {
                            std::string name_2 = atom_2_sel[isat]->name;
                            while (!name_2.empty() && name_2.back() == ' ') name_2.pop_back();
                            std::string pdb_atom_name_2 = coot::atom_id_mmdb_expand(name_2);

                            if (pdb_atom_name_2 == ltr.atom_id_2_4c()) {
                                for (int itat = 0; itat < n_atom_3; itat++) {
                                    std::string name_3 = atom_3_sel[itat]->name;
                                    while (!name_3.empty() && name_3.back() == ' ') name_3.pop_back();
                                    std::string pdb_atom_name_3 = coot::atom_id_mmdb_expand(name_3);

                                    if (pdb_atom_name_3 == ltr.atom_id_3_4c()) {
                                        for (int iffat = 0; iffat < n_atom_4; iffat++) {
                                            std::string name_4 = atom_4_sel[iffat]->name;
                                            while (!name_4.empty() && name_4.back() == ' ') name_4.pop_back();
                                            std::string pdb_atom_name_4 = coot::atom_id_mmdb_expand(name_4);

                                            if (pdb_atom_name_4 == ltr.atom_id_4_4c()) {

                                                // Skip dictionary mainchain torsions
                                                if (pdb_atom_name_1 == " N  " && pdb_atom_name_4 == " N  ") continue;
                                                if (pdb_atom_name_1 == " CA " && pdb_atom_name_4 == " CA ") continue;
                                                if (pdb_atom_name_1 == " C  " && pdb_atom_name_4 == " C  ") continue;

                                                // Check altloc compatibility
                                                char alt_conf_1 = atom_1_sel[ifat]->altloc;
                                                char alt_conf_2 = atom_2_sel[isat]->altloc;
                                                char alt_conf_3 = atom_3_sel[itat]->altloc;
                                                char alt_conf_4 = atom_4_sel[iffat]->altloc;

                                                // Normalization: gemmi uses ' ' for "no altloc"
                                                bool altloc_ok = ((alt_conf_1 == alt_conf_2) && (alt_conf_1 == alt_conf_3) && (alt_conf_1 == alt_conf_4));
                                                if (!altloc_ok) altloc_ok = ((alt_conf_1 == alt_conf_2) && (alt_conf_1 == alt_conf_3) && (alt_conf_4 == ' '));
                                                if (!altloc_ok) altloc_ok = ((alt_conf_1 == alt_conf_2) && (alt_conf_1 == alt_conf_4) && (alt_conf_3 == ' '));
                                                if (!altloc_ok) altloc_ok = ((alt_conf_1 == alt_conf_3) && (alt_conf_1 == alt_conf_4) && (alt_conf_2 == ' '));
                                                if (!altloc_ok) altloc_ok = ((alt_conf_2 == alt_conf_3) && (alt_conf_2 == alt_conf_4) && (alt_conf_1 == ' '));

                                                if (altloc_ok) {
                                                    n_torsions++;
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
        }
    }

    return n_torsions;
}

} // namespace coot