#pragma once

#include <string>
#include <gemmi/model.hpp>

namespace coot {

bool is_main_chain_p_gemmi(const gemmi::Atom& atom, const gemmi::Residue& residue) {
    const std::string& mol_atom_name = atom.name;
    if (mol_atom_name == "N" ||
        mol_atom_name == "C" ||
        mol_atom_name == "CA" ||
        mol_atom_name == "H" ||
        mol_atom_name == "HA" ||
        mol_atom_name == "OXT" ||
        mol_atom_name == "O") {
        return true;
    } else {
        const std::string& res_name = residue.name;
        if (res_name == "GLY") {
            if (mol_atom_name == "HA2" ||
                mol_atom_name == "HA3") {
                return true;
            }
        }
        return false;
    }
}

} // namespace coot