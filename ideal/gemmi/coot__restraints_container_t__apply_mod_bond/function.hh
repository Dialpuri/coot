#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>

#include "geometry/protein-geometry.hh"
#include "ideal/simple-restraint.hh"

namespace coot {

inline void apply_mod_bond_gemmi(std::vector<simple_restraint>& restraints_vec,
                                  const chem_mod_bond& mod_bond,
                                  gemmi::Residue& residue) {
    if (mod_bond.function == CHEM_MOD_FUNCTION_ADD) {
        int atom_count = static_cast<int>(residue.atoms.size());
        int idx1 = -1, idx2 = -1;
        for (int i = 0; i < atom_count; i++) {
            if (residue.atoms[static_cast<size_t>(i)].name == mod_bond.atom_id_1 && idx1 < 0)
                idx1 = i;
            if (residue.atoms[static_cast<size_t>(i)].name == mod_bond.atom_id_2 && idx2 < 0)
                idx2 = i;
        }
        if (idx1 >= 0 && idx2 >= 0) {
            simple_restraint restraint;
            restraint.restraint_type = BOND_RESTRAINT;
            restraint.atom_index_1 = idx1;
            restraint.atom_index_2 = idx2;
            restraint.target_value = mod_bond.new_value_dist;
            restraint.sigma = mod_bond.new_value_dist_esd;
            restraint.observed_value = 0.f;
            restraints_vec.push_back(restraint);
        }
    }
    if (mod_bond.function == CHEM_MOD_FUNCTION_CHANGE) {
        for (auto& r : restraints_vec) {
            if (r.restraint_type == BOND_RESTRAINT) {
                if (static_cast<size_t>(r.atom_index_1) < residue.atoms.size() &&
                    static_cast<size_t>(r.atom_index_2) < residue.atoms.size()) {
                    std::string name1 = residue.atoms[static_cast<size_t>(r.atom_index_1)].name;
                    std::string name2 = residue.atoms[static_cast<size_t>(r.atom_index_2)].name;
                    if ((name1 == mod_bond.atom_id_1 && name2 == mod_bond.atom_id_2) ||
                        (name2 == mod_bond.atom_id_1 && name1 == mod_bond.atom_id_2)) {
                        r.target_value = mod_bond.new_value_dist;
                        r.sigma = mod_bond.new_value_dist_esd;
                    }
                }
            }
        }
    }
    if (mod_bond.function == CHEM_MOD_FUNCTION_DELETE) {
        bool atom1_in_res = false;
        bool atom2_in_res = false;
        for (auto& atom : residue.atoms) {
            if (atom.name == mod_bond.atom_id_1) atom1_in_res = true;
            if (atom.name == mod_bond.atom_id_2) atom2_in_res = true;
        }

        if (!atom1_in_res || !atom2_in_res) return;

        int i = 0;
        while (i < static_cast<int>(restraints_vec.size())) {
            if (restraints_vec[i].restraint_type == BOND_RESTRAINT) {
                if (static_cast<size_t>(restraints_vec[i].atom_index_1) < residue.atoms.size() &&
                    static_cast<size_t>(restraints_vec[i].atom_index_2) < residue.atoms.size()) {
                    std::string name1 = residue.atoms[static_cast<size_t>(restraints_vec[i].atom_index_1)].name;
                    std::string name2 = residue.atoms[static_cast<size_t>(restraints_vec[i].atom_index_2)].name;
                    if ((name1 == mod_bond.atom_id_1 && name2 == mod_bond.atom_id_2) ||
                        (name2 == mod_bond.atom_id_1 && name1 == mod_bond.atom_id_2)) {
                        restraints_vec[i] = restraints_vec.back();
                        restraints_vec.pop_back();
                        continue;
                    }
                }
            }
            i++;
        }
    }
}

} // namespace coot