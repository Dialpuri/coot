#include "function.hh"
#include <vector>
#include <set>
#include <string>
#include <algorithm>

namespace coot {

int add_angles_gemmi(
    std::vector<simple_restraint>& restraints_vec,
    std::vector<std::set<int>>& bonded_atom_indices,
    std::set<int>& fixed_atom_indices,
    int idr,
    const std::vector<gemmi::Atom*>& res_atoms,
    gemmi::Residue* sel_res,
    const coot::protein_geometry& geom) {

    int n_angle_restr = 0;

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

    int i_no_res_atoms = static_cast<int>(res_atoms.size());

    // Build string_atom_names array
    std::vector<std::string> string_atom_names(i_no_res_atoms);
    for (int iat = 0; iat < i_no_res_atoms; iat++) {
        string_atom_names[iat] = pad_to_4(res_atoms[iat]->name);
    }

    for (unsigned int ib = 0; ib < dict.angle_restraint.size(); ib++) {
        for (int iat = 0; iat < i_no_res_atoms; iat++) {
            const std::string& pdb_atom_name1 = string_atom_names[iat];

            if (pdb_atom_name1 == dict.angle_restraint[ib].atom_id_1_4c()) {
                for (int iat2 = 0; iat2 < i_no_res_atoms; iat2++) {
                    const std::string& pdb_atom_name2 = string_atom_names[iat2];

                    if (pdb_atom_name2 == dict.angle_restraint[ib].atom_id_2_4c()) {
                        for (int iat3 = 0; iat3 < i_no_res_atoms; iat3++) {
                            const std::string& pdb_atom_name3 = string_atom_names[iat3];

                            if (pdb_atom_name3 == dict.angle_restraint[ib].atom_id_3_4c()) {
                                // Check alt locs - in gemmi, altloc is a char, '\0' means no altloc
                                char alt_1 = res_atoms[iat]->altloc;
                                char alt_2 = res_atoms[iat2]->altloc;
                                char alt_3 = res_atoms[iat3]->altloc;

                                // In MMDB, empty string "" means no altloc.
                                // In gemmi, '\0' means no altloc.
                                // The original logic: all three match, OR any is "" and the other two match
                                if (((alt_1 == alt_2) && (alt_1 == alt_3)) ||
                                    ((alt_1 == '\0') && (alt_2 == alt_3)) ||
                                    ((alt_2 == '\0') && (alt_1 == alt_3)) ||
                                    ((alt_3 == '\0') && (alt_1 == alt_2))) {

                                    int index1 = iat;
                                    int index2 = iat2;
                                    int index3 = iat3;

                                    // Ensure bonded_atom_indices is sized properly
                                    if (index1 >= static_cast<int>(bonded_atom_indices.size()))
                                        bonded_atom_indices.resize(index1 + 1);
                                    if (index3 >= static_cast<int>(bonded_atom_indices.size()))
                                        bonded_atom_indices.resize(index3 + 1);

                                    bonded_atom_indices[index1].insert(index3);
                                    bonded_atom_indices[index3].insert(index1);

                                    // make_fixed_flags: for angle restraints, 
                                    // fixed atoms are those in fixed_atom_indices
                                    // that are NOT the middle (flanking) atom
                                    std::vector<bool> fixed_flags(3, false);
                                    if ((index1 != index2) && fixed_atom_indices.count(index1))
                                        fixed_flags[0] = true;
                                    if ((index2 != index1) && fixed_atom_indices.count(index2))
                                        fixed_flags[1] = true;
                                    if ((index3 != index2) && fixed_atom_indices.count(index3))
                                        fixed_flags[2] = true;

                                    // is_single_H_atom_angle_restraint: true if exactly one of
                                    // atom1 or atom3 (the outer atoms) is hydrogen
                                    bool is_H1 = res_atoms[iat]->element.is_hydrogen();
                                    bool is_H3 = res_atoms[iat3]->element.is_hydrogen();
                                    bool is_single_H_atom_angle_restraint = (is_H1 ^ is_H3);

                                    restraints_vec.push_back(simple_restraint(
                                        ANGLE_RESTRAINT,
                                        index1,
                                        index2,
                                        index3,
                                        fixed_flags,
                                        dict.angle_restraint[ib].angle(),
                                        dict.angle_restraint[ib].esd(),
                                        is_single_H_atom_angle_restraint
                                    ));
                                    n_angle_restr++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return n_angle_restr;
}

} // namespace coot