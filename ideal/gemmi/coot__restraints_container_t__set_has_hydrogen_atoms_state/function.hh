#pragma once

#include <gemmi/model.hpp>
#include <string>

namespace coot {

// Helper function: check if an atom is a hydrogen (or deuterium)
// This is the gemmi version of the original is_hydrogen function
inline bool
is_hydrogen(const gemmi::Atom& atom) {
    std::string ele = atom.element.name();
    return (ele == "H" || ele == " D" || ele == "D");
}

// Port of set_has_hydrogen_atoms_state to gemmi-style API
// This function checks if any atom in the given atoms vector is a hydrogen
// and updates the model_has_hydrogen_atoms flag accordingly
void set_has_hydrogen_atoms_state(const std::vector<gemmi::Atom>& atoms, bool& model_has_hydrogen_atoms);

} // namespace coot