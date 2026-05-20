#pragma once

#include <gemmi/model.hpp>

namespace coot {

// Port of MMDB's is_hydrogen_atom to gemmi
// Checks if an atom is a hydrogen (H) or deuterium (D) atom
inline bool is_hydrogen_atom_gemmi(const gemmi::Atom& atom) {
    return atom.is_hydrogen();
}

}