#pragma once
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>

namespace coot {

// Remove TER atoms from a residue.
// In gemmi, TER atoms are identified by having "TER" in their name.
void remove_TER_internal_gemmi(gemmi::Structure& st, gemmi::Residue* residue_p) {
    bool deleted = false;
    if (residue_p) {
        for (auto it = residue_p->atoms.begin(); it != residue_p->atoms.end();) {
            // Check if atom name contains "TER" (like "TER ")
            if (it->name.find("TER") != std::string::npos) {
                it = residue_p->atoms.erase(it);
                deleted = true;
            } else {
                ++it;
            }
        }
    }
}

} // namespace coot