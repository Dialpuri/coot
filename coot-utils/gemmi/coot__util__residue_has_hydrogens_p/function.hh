#pragma once
#include <gemmi/model.hpp>

namespace coot { namespace util {

// Port of coot::util::residue_has_hydrogens_p from MMDB to gemmi.
// Returns true if the given residue contains any hydrogen or deuterium atoms.
inline bool
residue_has_hydrogens_p_gemmi(const gemmi::Residue& residue) {
    for (const gemmi::Atom& atom : residue.atoms) {
        if (atom.element.is_hydrogen()) {
            return true;
        }
    }
    return false;
}

}} // namespace coot::util