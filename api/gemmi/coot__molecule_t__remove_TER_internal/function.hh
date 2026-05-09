#pragma once
#include <gemmi/model.hpp>

namespace coot { namespace molecule_t {

void remove_TER_internal_gemmi(gemmi::Residue& residue) {
    // Remove atoms that are TER atoms from the residue
    for (int i = static_cast<int>(residue.atoms.size()) - 1; i >= 0; --i) {
        // In PDB format, TER atoms have names containing "TER"
        if (residue.atoms[i].name.find("TER") != std::string::npos) {
            residue.atoms.erase(residue.atoms.begin() + i);
        }
    }
}

}}