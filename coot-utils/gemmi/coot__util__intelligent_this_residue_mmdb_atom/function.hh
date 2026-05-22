#pragma once
#include <gemmi/model.hpp>
#include <string>

namespace coot { namespace util {

inline const gemmi::Atom* intelligent_this_residue_mmdb_atom_gemmi(const gemmi::Residue& residue) {
    for (const auto& atom : residue.atoms) {
        if (atom.name == "CA" || atom.name == " C1'" || atom.name == " CA ") {
            return &atom;
        }
    }
    
    if (!residue.atoms.empty()) {
        return &residue.atoms[0];
    }
    
    return nullptr;
}

}} // namespace coot::util