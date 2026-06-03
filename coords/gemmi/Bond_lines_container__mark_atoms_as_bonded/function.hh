#pragma once
#include <gemmi/model.hpp>
#include <map>
#include <string>

// gemmi replacement for graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND
namespace graphical_bonds_container {
    constexpr int BONDED_WITH_STANDARD_ATOM_BOND = 1;
}

// gemmi port of Bond_lines_container::mark_atoms_as_bonded
// Uses std::map<gemmi::Atom*, int> instead of MMDB UDData on atoms.
//
// Semantics: when have_udd_atoms==true and done_bond_udd_handle==false,
// mark each atom (unless its element is S/Se/Cl/Br/P) in the bonded_map
// with value BONDED_WITH_STANDARD_ATOM_BOND.
inline void mark_atoms_as_bonded_gemmi(
    gemmi::Atom *atom_p_1, gemmi::Atom *atom_p_2,
    bool have_udd_atoms,
    std::map<gemmi::Atom*, int>& bonded_map,
    bool done_bond_udd_handle)
{
    auto is_excluded = [](const std::string& elem) {
        return elem == "S"  || elem == "Se" ||
               elem == "Cl" || elem == "Br" || elem == "P";
    };

    if (have_udd_atoms && !done_bond_udd_handle) {
        if (!is_excluded(atom_p_1->element.name())) {
            bonded_map[atom_p_1] = graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND;
        }
        if (!is_excluded(atom_p_2->element.name())) {
            bonded_map[atom_p_2] = graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND;
        }
    }
}