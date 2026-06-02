#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <gemmi/model.hpp>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__torsion_general__atom_index/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__torsion_general__get_contact_indices/gemmi/function.hh"

namespace coot {
namespace torsion_general {

struct torsion_general_result {
    bool setup_correctly;
    gemmi::Residue* residue_p;
    std::vector<atom_spec_t> user_defined_torsion_atoms;
    std::vector<std::vector<int>> contact_indices;
    std::vector<int> clicked_atom_indices;

    torsion_general_result() : setup_correctly(false), residue_p(nullptr) {}
};

inline torsion_general_result torsion_general_gemmi(
    gemmi::Residue* res,
    const std::vector<atom_spec_t>& user_defined_torsion_atoms_in) {

    torsion_general_result result;
    result.setup_correctly = false;
    result.residue_p = res;
    result.user_defined_torsion_atoms = user_defined_torsion_atoms_in;

    if (!res) return result;

    if (user_defined_torsion_atoms_in.size() < 4) return result;

    int index_0 = atom_index_gemmi(user_defined_torsion_atoms_in[0], *res);
    int index_1 = atom_index_gemmi(user_defined_torsion_atoms_in[1], *res);
    int index_2 = atom_index_gemmi(user_defined_torsion_atoms_in[2], *res);
    int index_3 = atom_index_gemmi(user_defined_torsion_atoms_in[3], *res);

    const int FAIL = -1;

    if (index_0 != FAIL) {
        if (index_1 != FAIL) {
            if (index_2 != FAIL) {
                if (index_3 != FAIL) {
                    result.contact_indices = get_contact_indices_gemmi(*res);
                    result.clicked_atom_indices.clear();
                    result.clicked_atom_indices.push_back(index_0);
                    result.clicked_atom_indices.push_back(index_1);
                    result.clicked_atom_indices.push_back(index_2);
                    result.clicked_atom_indices.push_back(index_3);
                    result.setup_correctly = true;
                } else {
                    std::cout << "ERROR:: failed to find " << user_defined_torsion_atoms_in[3].atom_name << std::endl;
                }
            } else {
                std::cout << "ERROR:: failed to find " << user_defined_torsion_atoms_in[2].atom_name << std::endl;
            }
        } else {
            std::cout << "ERROR:: failed to find " << user_defined_torsion_atoms_in[1].atom_name << std::endl;
        }
    } else {
        std::cout << "ERROR:: failed to find " << user_defined_torsion_atoms_in[0].atom_name << std::endl;
    }

    return result;
}

} // namespace torsion_general
} // namespace coot