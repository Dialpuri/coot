#pragma once

#include <string>
#include <gemmi/model.hpp>

namespace coot {

// String overload: checks space-padded 4-character PDB atom names.
// Note: includes " CB " (different from the atom overload).
inline bool is_main_chain_p_gemmi(const std::string &mol_atom_name) {
    if (mol_atom_name == " N  " ||
        mol_atom_name == " C  " ||
        mol_atom_name == " H  " ||
        mol_atom_name == " CA " ||
        mol_atom_name == " CB " ||
        mol_atom_name == " HA " || // CA hydrogen
        mol_atom_name == " O  ") {
        return true;
    }
    return false;
}

// Atom overload: checks unpadded atom name (gemmi convention) and
// handles GLY HA2/HA3 as main-chain. Excludes CB (sidechain).
// Requires residue reference because gemmi::Atom has no parent pointer.
inline bool is_main_chain_p_gemmi(const gemmi::Atom &atom, const gemmi::Residue &residue) {
    const std::string &mol_atom_name = atom.name;
    if (mol_atom_name == "N" ||
        mol_atom_name == "C" ||
        mol_atom_name == "CA" ||
        mol_atom_name == "H" ||
        mol_atom_name == "HA" ||
        mol_atom_name == "OXT" ||
        mol_atom_name == "O") {
        return true;
    }
    if (residue.name == "GLY") {
        if (mol_atom_name == "HA2" ||
            mol_atom_name == "HA3") {
            return true;
        }
    }
    return false;
}

} // namespace coot
