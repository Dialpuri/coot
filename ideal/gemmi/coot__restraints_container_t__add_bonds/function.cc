#include "function.hh"
#include <vector>
#include <string>
#include <algorithm>

namespace coot {

int add_bonds_gemmi(
    std::vector<simple_restraint>& restraints_vec,
    std::vector<std::set<int>>& bonded_atom_indices,
    std::map<gemmi::Atom*, hb_t>& H_atom_parent_energy_type_atom_map,
    int idr,
    const std::vector<gemmi::Atom*>& res_atoms,
    gemmi::Residue* sel_res,
    const coot::protein_geometry& geom) {

    int n_bond_restr = 0;

    const dictionary_residue_restraints_t& dict = geom[idr].second;

    // Helper: pad atom name to 4 characters matching dictionary convention
    // All names get a leading space then right-padded to 4.
    // Examples: "N" -> " N  ", "CA" -> " CA ", "HB3" -> " HB3", "OXT" -> " OXT"
    auto pad_to_4 = [](const std::string& name) -> std::string {
        std::string result = " " + name;
        if (result.size() < 4)
            result.resize(4, ' ');
        return result;
    };

    for (unsigned int ib = 0; ib < dict.bond_restraint.size(); ib++) {
        for (size_t iat = 0; iat < res_atoms.size(); iat++) {
            std::string pdb_atom_name1 = pad_to_4(res_atoms[iat]->name);

            if (pdb_atom_name1 == dict.bond_restraint[ib].atom_id_1_4c()) {
                for (size_t iat2 = 0; iat2 < res_atoms.size(); iat2++) {
                    std::string pdb_atom_name2 = pad_to_4(res_atoms[iat2]->name);

                    if (pdb_atom_name2 == dict.bond_restraint[ib].atom_id_2_4c()) {
                        // Check alt confs - in gemmi, altloc is a char, '\0' means no altloc
                        char alt_1 = res_atoms[iat]->altloc;
                        char alt_2 = res_atoms[iat2]->altloc;
                        if (alt_1 == '\0' || alt_2 == '\0' || alt_1 == alt_2) {
                            // Atom indices
                            int index1 = static_cast<int>(iat);
                            int index2 = static_cast<int>(iat2);

                            // Ensure bonded_atom_indices is sized properly
                            if (index1 >= static_cast<int>(bonded_atom_indices.size()))
                                bonded_atom_indices.resize(index1 + 1);
                            if (index2 >= static_cast<int>(bonded_atom_indices.size()))
                                bonded_atom_indices.resize(index2 + 1);

                            bonded_atom_indices[index1].insert(index2);
                            bonded_atom_indices[index2].insert(index1);

                            // Fixed flags - both false for monomer bonds
                            std::vector<bool> fixed_flags = {false, false};

                            try {
                                restraints_vec.push_back(simple_restraint(
                                    BOND_RESTRAINT,
                                    index1,
                                    index2,
                                    fixed_flags,
                                    dict.bond_restraint[ib].value_dist(),
                                    dict.bond_restraint[ib].value_esd(),
                                    1.2
                                ));
                                n_bond_restr++;

                                // H-atom parent energy type tracking
                                if (res_atoms[iat]->element.is_hydrogen()) {
                                    gemmi::Atom* H_at = res_atoms[iat];
                                    std::string parent_atom_name(res_atoms[iat2]->name);
                                    std::string te = dict.type_energy(parent_atom_name);
                                    hb_t hbt = geom.get_h_bond_type(te);
                                    H_atom_parent_energy_type_atom_map[H_at] = hbt;

                                    // HIS kludge
                                    if (dict.residue_info.comp_id == "HIS") {
                                        if (parent_atom_name == "ND1 ")
                                            H_atom_parent_energy_type_atom_map[H_at] = HB_BOTH;
                                        if (parent_atom_name == "NE2 ")
                                            H_atom_parent_energy_type_atom_map[H_at] = HB_BOTH;
                                    }
                                    // TRP kludge
                                    if (dict.residue_info.comp_id == "TRP") {
                                        if (parent_atom_name == "NE1 ")
                                            H_atom_parent_energy_type_atom_map[H_at] = HB_BOTH;
                                    }
                                }

                                if (res_atoms[iat2]->element.is_hydrogen()) {
                                    gemmi::Atom* H_at = res_atoms[iat2];
                                    std::string parent_atom_name(res_atoms[iat]->name);
                                    std::string te = dict.type_energy(parent_atom_name);
                                    hb_t hbt = geom.get_h_bond_type(te);
                                    H_atom_parent_energy_type_atom_map[H_at] = hbt;

                                    // HIS kludge
                                    if (dict.residue_info.comp_id == "HIS") {
                                        if (parent_atom_name == "ND1 ")
                                            H_atom_parent_energy_type_atom_map[H_at] = HB_BOTH;
                                        if (parent_atom_name == "NE2 ")
                                            H_atom_parent_energy_type_atom_map[H_at] = HB_BOTH;
                                    }
                                    // TRP kludge
                                    if (dict.residue_info.comp_id == "TRP") {
                                        if (parent_atom_name == "NE1 ")
                                            H_atom_parent_energy_type_atom_map[H_at] = HB_BOTH;
                                    }
                                }
                            } catch (const std::runtime_error& rte) {
                                // Silently ignore missing target geometry
                                std::cout << "trapped a runtime_error on adding bond restraint: "
                                          << rte.what() << std::endl;
                            }
                        }
                    }
                }
            }
        }
    }

    return n_bond_restr;
}

} // namespace coot