#include "function.hh"

namespace coot {

void mod_bond_change_gemmi(std::vector<simple_restraint>& restraints_vec,
                           const std::vector<std::string>& atom_names,
                           const std::vector<gemmi::Residue*>& atom_residues,
                           const chem_mod_bond& mod_bond,
                           gemmi::Residue* target_residue) {
    for (unsigned int i = 0; i < restraints_vec.size(); i++) {
        simple_restraint& rest = restraints_vec[i];
        if (rest.restraint_type == BOND_RESTRAINT) {
            if (atom_residues[rest.atom_index_1] == target_residue) {
                if (atom_residues[rest.atom_index_2] == target_residue) {
                    std::string name_1 = atom_names[rest.atom_index_1];
                    std::string name_2 = atom_names[rest.atom_index_2];
                    if (name_1 == mod_bond.atom_id_1) {
                        if (name_2 == mod_bond.atom_id_2) {
                            rest.target_value = mod_bond.new_value_dist;
                            rest.sigma = mod_bond.new_value_dist_esd;
                        }
                    }
                }
            }
        }
    }
}

}  // namespace coot