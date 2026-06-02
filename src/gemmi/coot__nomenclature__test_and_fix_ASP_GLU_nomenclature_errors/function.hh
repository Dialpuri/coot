#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <iostream>
#include <utility>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__atom_index_quad__torsion__10c420/gemmi/function.hh"

namespace coot {
namespace nomenclature {

int test_and_fix_ASP_GLU_nomenclature_errors_gemmi(
    gemmi::Residue &residue, bool apply_swap_if_found) {

    int iswapped = 0;
    int n_atoms = int(residue.atoms.size());

    std::string residue_name = residue.name;

    if (residue_name == "ASP" || residue_name == "GLU") {

        // gemmi atom names are unpadded (e.g. "OD1" not " OD1")
        std::string test_atom_name = "OD1";
        if (residue_name == "GLU")
            test_atom_name = "OE1";

        std::vector<char> alt_conf_list;
        // first get the altconfs in the residue:
        for (int i = 0; i < n_atoms; i++) {
            if (residue.atoms[i].name == test_atom_name) {
                alt_conf_list.push_back(residue.atoms[i].altloc);
            }
        }

        for (unsigned int ialtconf = 0; ialtconf < alt_conf_list.size(); ialtconf++) {

            coot::atom_index_quad quad(-1, -1, -1, -1);

            for (int i = 0; i < n_atoms; i++) {
                if (residue.atoms[i].altloc == alt_conf_list[ialtconf]) {

                    if (residue_name == "ASP") {
                        if (residue.atoms[i].name == "CA") quad.index1 = i;
                        if (residue.atoms[i].name == "CB") quad.index2 = i;
                        if (residue.atoms[i].name == "CG") quad.index3 = i;
                        if (residue.atoms[i].name == "OD1") quad.index4 = i;
                    }

                    if (residue_name == "GLU") {
                        if (residue.atoms[i].name == "CB") quad.index1 = i;
                        if (residue.atoms[i].name == "CG") quad.index2 = i;
                        if (residue.atoms[i].name == "CD") quad.index3 = i;
                        if (residue.atoms[i].name == "OE1") quad.index4 = i;
                    }
                }
            }

            try {
                double torsion = quad.torsion_gemmi(residue);
                bool torsion_is_good = false;
                if (torsion >= -90.0) {
                    if (torsion <= 90.0) {
                        torsion_is_good = true;
                    }
                }
                if (!torsion_is_good) {

                    // SWAP names (unpadded gemmi names)
                    std::string swap_name_1 = "OD1";
                    std::string swap_name_2 = "OD2";
                    if (residue_name == "GLU") {
                        swap_name_1 = "OE1";
                        swap_name_2 = "OE2";
                    }

                    int at_1_idx = -1;
                    int at_2_idx = -1;
                    for (int i = 0; i < n_atoms; i++) {
                        if (residue.atoms[i].altloc == alt_conf_list[ialtconf]) {
                            if (residue.atoms[i].name == swap_name_1) at_1_idx = i;
                            if (residue.atoms[i].name == swap_name_2) at_2_idx = i;
                        }
                    }
                    if (at_1_idx >= 0 && at_2_idx >= 0) {
                        if (apply_swap_if_found) {
                            std::swap(residue.atoms[at_1_idx].name,
                                      residue.atoms[at_2_idx].name);
                        }
                        iswapped = 1;
                    }
                }
            } catch (const std::runtime_error &rte) {
                std::cout << "WARNING:: missing atoms " << rte.what() << std::endl;
            }
        }
    }

    return iswapped;
}

} // namespace nomenclature
} // namespace coot