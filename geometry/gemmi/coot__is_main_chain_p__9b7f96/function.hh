#pragma once
#include <gemmi/model.hpp>
#include <string>

namespace coot {

inline bool is_main_chain_p_gemmi(gemmi::CRA cra) {
    if (!cra.atom || !cra.residue || !cra.chain) return false;

    const gemmi::Atom* at = cra.atom;
    std::string mol_atom_name = at->name;
    mol_atom_name.erase(mol_atom_name.find_last_not_of(" \t\r\n") + 1);
    mol_atom_name.erase(0, mol_atom_name.find_first_not_of(" \t\r\n"));

    if (mol_atom_name == "N" ||
        mol_atom_name == "C" ||
        mol_atom_name == "CA" ||
        mol_atom_name == "H" ||
        mol_atom_name == "HA" ||
        mol_atom_name == "OXT" ||
        mol_atom_name == "O") {
        return true;
    }

    std::string res_name = cra.residue->name;
    res_name.erase(res_name.find_last_not_of(" \t\r\n") + 1);
    res_name.erase(0, res_name.find_first_not_of(" \t\r\n"));

    if (res_name == "GLY") {
        if (mol_atom_name == "HA2" || mol_atom_name == "HA3") {
            return true;
        }

        // N-terminus approximation: first residue in its chain
        if (cra.residue == &cra.chain->residues.front()) {
            if (mol_atom_name == "H1") return true;
            if (mol_atom_name == "H2") return true;
            if (mol_atom_name == "H3") return true;
        }
    }

    return false;
}

} // namespace coot
