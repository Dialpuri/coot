#pragma once
#include <string>
#include <gemmi/model.hpp>

namespace coot {
namespace ideal_rna {

void fix_up_residue_and_atom_names_gemmi(gemmi::Residue* residue, bool is_dna_flag) {
    std::string res_name = residue->name;
    std::string new_name = res_name; // residue_name_old_to_new(res_name, is_dna_flag); no longer needed

    residue->name = new_name;

    for (auto& atom : residue->atoms) {
        std::string atom_name = atom.name;
        if (atom_name.length() > 3) {
            if (atom_name[3] == '*') {
                atom_name[3] = '\'';
                atom.name = atom_name;
            }

            if (atom_name == " O1P") atom.name = " OP1";
            if (atom_name == " O2P") atom.name = " OP2";
        }
    }

    // fix the atom name C5M->C7 on a T in DNA [Grr, !@#$!@#$% PDB...]
    if (new_name == "DT") {
        for (auto& atom : residue->atoms) {
            std::string atom_name = atom.name;
            if (atom_name == " C5M") {
                atom.name = " C7 ";
            }
        }
    }
}

} // namespace ideal_rna
} // namespace coot