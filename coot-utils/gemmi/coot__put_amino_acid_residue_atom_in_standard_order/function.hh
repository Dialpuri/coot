#pragma once
#include <gemmi/model.hpp>
#include <vector>
#include <string>

namespace coot {

void put_amino_acid_residue_atom_in_standard_order_gemmi(gemmi::Residue& residue) {
    std::vector<gemmi::Atom*> N_ats;
    std::vector<gemmi::Atom*> other_ats;
    
    for (auto& atom : residue.atoms) {
        if (atom.name == " N  ") { // PDBv3 format
            N_ats.push_back(&atom);
        } else {
            other_ats.push_back(&atom);
        }
    }
    
    // Reorder: N atoms first, then others
    size_t idx = 0;
    for (auto* atom : N_ats) {
        residue.atoms[idx] = *atom;
        idx++;
    }
    for (auto* atom : other_ats) {
        residue.atoms[idx] = *atom;
        idx++;
    }
}

} // namespace coot