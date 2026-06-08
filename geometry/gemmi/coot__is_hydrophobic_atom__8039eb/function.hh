#pragma once
#include <string>

namespace coot {

inline bool is_hydrophobic_atom_gemmi(const std::string &residue_name, const std::string &atom_name) {
    bool status = false;

    if (atom_name.size() > 1 && atom_name[1] == 'N') return false;
    if (atom_name.size() > 1 && atom_name[1] == 'O') return false;

    if (residue_name == "GLY")
        status = true;
    else if (residue_name == "ALA")
        status = true;
    else if (residue_name == "VAL")
        status = true;
    else if (residue_name == "LEU")
        status = true;
    else if (residue_name == "ILE")
        status = true;
    else if (residue_name == "PRO")
        status = true;
    else if (residue_name == "PHE")
        status = true;
    else if (residue_name == "MET")
        status = true;
    else if (residue_name == "TRP")
        status = true;
    else if (residue_name == "CYS")
        status = true;
    else if (residue_name == "TYR")
        status = true;

    return status;
}

} // namespace coot