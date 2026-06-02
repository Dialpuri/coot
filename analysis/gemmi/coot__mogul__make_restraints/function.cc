#include "function.hh"

#include <vector>
#include <string>

namespace coot {

coot::dictionary_residue_restraints_t
make_restraints_gemmi(
    const gemmi::Residue* residue_p,
    const std::string& comp_id,
    const std::vector<coot::mogul_item>& items,
    int imol,
    const coot::protein_geometry& geom
) {
    coot::dictionary_residue_restraints_t r(comp_id, -1);

    if (residue_p) {
        // Get atom names from gemmi::Residue
        std::vector<std::string> atom_names;
        for (const auto& atom : residue_p->atoms) {
            atom_names.push_back(atom.name);
        }

        int n_residue_atoms = static_cast<int>(atom_names.size());

        std::pair<bool, dictionary_residue_restraints_t> current_restraints =
            geom.get_monomer_restraints(comp_id, imol);

        for (size_t i = 0; i < items.size(); i++) {
            if (items[i].type == coot::mogul_item::BOND) {
                int idx_1 = items[i].idx_1 - 1;
                int idx_2 = items[i].idx_2 - 1;
                if (idx_1 >= 0 && idx_1 < n_residue_atoms) {
                    if (idx_2 >= 0 && idx_2 < n_residue_atoms) {
                        std::string name_1 = atom_names[idx_1];
                        std::string name_2 = atom_names[idx_2];
                        std::string type;
                        if (current_restraints.first)
                            type = current_restraints.second.get_bond_type(name_1, name_2);
                        float dist = items[i].median;
                        float esd  = items[i].std_dev;
                        dict_bond_restraint_t rest(name_1, name_2, type, dist, esd, 0.0, 0.0, false);
                        r.bond_restraint.push_back(rest);
                    }
                }
            }

            if (items[i].type == coot::mogul_item::ANGLE) {
                int idx_1 = items[i].idx_1 - 1;
                int idx_2 = items[i].idx_2 - 1;
                int idx_3 = items[i].idx_3 - 1;
                if (idx_1 >= 0 && idx_1 < n_residue_atoms) {
                    if (idx_2 >= 0 && idx_2 < n_residue_atoms) {
                        if (idx_3 >= 0 && idx_3 < n_residue_atoms) {
                            std::string name_1 = atom_names[idx_1];
                            std::string name_2 = atom_names[idx_2];
                            std::string name_3 = atom_names[idx_3];
                            float angle = items[i].median;
                            float esd   = items[i].std_dev;
                            dict_angle_restraint_t rest(name_1, name_2, name_3, angle, esd);
                            r.angle_restraint.push_back(rest);
                        }
                    }
                }
            }
        }
    }
    return r;
}

} // namespace coot