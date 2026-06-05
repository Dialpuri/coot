#include "function.hh"
#include <iostream>
#include <algorithm>

namespace coot {

int add_link_bond_gemmi(
    std::vector<simple_restraint>& restraints_vec,
    std::vector<std::set<int>>& bonded_atom_indices,
    std::string link_type,
    gemmi::Residue const *first,
    gemmi::Residue const *second,
    short int is_fixed_first,
    short int is_fixed_second,
    const coot::protein_geometry &geom,
    int first_residue_index,
    int second_residue_index) {

    bool debug = false;

    int n_first_res_atoms  = static_cast<int>(first->atoms.size());
    int n_second_res_atoms = static_cast<int>(second->atoms.size());

    std::vector<bool> fixed_atom_flags(2); // 2 atoms in this restraint

    short int found_link_type = 0;
    int index1, index2;

    if (n_first_res_atoms <= 0) {
        std::cout << "no atoms in first residue!? " << std::endl;
    }
    if (n_second_res_atoms <= 0) {
        std::cout << "no atoms in second residue!? " << std::endl;
    }

    if (debug) {
        std::cout << "INFO:: geom.link_size() is " << geom.link_size() << std::endl;
        std::cout << "first residue:\n";
        for (int i = 0; i < n_first_res_atoms; i++)
            std::cout << "    " << first->atoms[i].name
                      << " " << first->seqid.num.value << "\n";
        std::cout << "second residue:\n";
        for (int i = 0; i < n_second_res_atoms; i++)
            std::cout << "    " << second->atoms[i].name
                      << " " << second->seqid.num.value << "\n";
    }

    int nbond = 0;

    // Helper: pad atom name to 4 characters matching dictionary convention.
    // All names get a leading space then right-padded to 4.
    auto pad_to_4 = [](const std::string& name) -> std::string {
        std::string result = " " + name;
        if (result.size() < 4)
            result.resize(4, ' ');
        return result;
    };

    // Helper: check if an atom index is already involved in a "fixed" restraint
    // (i.e. appears in any existing simple_restraint). We scan restraints_vec.
    auto is_atom_fixed = [&](int idx) -> bool {
        for (const auto& restr : restraints_vec) {
            if (restr.atom_index_1 == idx) return true;
            if (restr.atom_index_2 == idx) return true;
        }
        return false;
    };

    for (int i = 0; i < geom.link_size(); i++) {
        if (geom.link(i).link_id == link_type) { // typically "TRANS"
            found_link_type = 1;
            for (unsigned int j = 0; j < geom.link(i).link_bond_restraint.size(); j++) {
                if (geom.link(i).link_bond_restraint[j].atom_1_comp_id == 1 &&
                    geom.link(i).link_bond_restraint[j].atom_2_comp_id == 2) {
                    // as expected
                } else {
                    std::cout << "PROGRAMMER ERROR (shortsighted fool)" << std::endl;
                    std::cout << "bad things will now happen..." << std::endl;
                }

                for (int ifat = 0; ifat < n_first_res_atoms; ifat++) {
                    std::string pdb_atom_name_1 = pad_to_4(first->atoms[ifat].name);

                    if (pdb_atom_name_1 == geom.link(i).link_bond_restraint[j].atom_id_1_4c()) {
                        for (int isat = 0; isat < n_second_res_atoms; isat++) {
                            std::string pdb_atom_name_2 = pad_to_4(second->atoms[isat].name);

                            if (pdb_atom_name_2 == geom.link(i).link_bond_restraint[j].atom_id_2_4c()) {

                                if (debug)
                                    std::cout << "DEBUG::  adding " << link_type
                                              << " bond for "
                                              << first->seqid.num.value
                                              << " -> " << second->seqid.num.value
                                              << " atoms "
                                              << first->atoms[ifat].name
                                              << " to "
                                              << second->atoms[isat].name
                                              << std::endl;

                                // Check that alt confs match
                                // In gemmi, altloc is a char; '\0' means no alternate location
                                char alt_conf_1  = first->atoms[ifat].altloc;
                                char alt_conf_2  = second->atoms[isat].altloc;
                                if ((alt_conf_1 == alt_conf_2) ||
                                    (alt_conf_1 == '\0') ||
                                    (alt_conf_2 == '\0')) {

                                    // Compute global atom indices
                                    index1 = first_residue_index * 100 + ifat;
                                    index2 = second_residue_index * 100 + isat;

                                    // Resize bonded_atom_indices if needed
                                    if (index1 >= static_cast<int>(bonded_atom_indices.size()))
                                        bonded_atom_indices.resize(index1 + 1);
                                    if (index2 >= static_cast<int>(bonded_atom_indices.size()))
                                        bonded_atom_indices.resize(index2 + 1);

                                    bonded_atom_indices[index1].insert(index2);
                                    bonded_atom_indices[index2].insert(index1);

                                    fixed_atom_flags[0] = is_fixed_first;
                                    fixed_atom_flags[1] = is_fixed_second;

                                    // Check existing fixed flags from other restraints
                                    // (equivalent to restraints_container_t::make_fixed_flags)
                                    if (is_atom_fixed(index1)) fixed_atom_flags[0] = 1;
                                    if (is_atom_fixed(index2)) fixed_atom_flags[1] = 1;

                                    restraints_vec.push_back(simple_restraint(
                                        BOND_RESTRAINT,
                                        index1,
                                        index2,
                                        fixed_atom_flags,
                                        geom.link(i).link_bond_restraint[j].dist(),
                                        geom.link(i).link_bond_restraint[j].esd(),
                                        1.2 // junk value
                                    ));
                                    nbond++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if (found_link_type == 0)
        std::cout << "link type \"" << link_type << "\" not found in dictionary!!\n";

    return nbond;
}

} // namespace coot