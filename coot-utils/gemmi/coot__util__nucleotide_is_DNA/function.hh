#pragma once
#include <string>
#include <gemmi/model.hpp>

namespace coot {
namespace util {

// Original MMDB function: checks if residue is DNA (no O2' atom present)
// Returns true if residue has no O2' or O2* atom (DNA logic)
// Returns false if residue has O2' or O2* atom (RNA logic)
inline bool nucleotide_is_DNA_gemmi(const gemmi::Residue* residue) {
    if (!residue) {
        // nullptr residue: no atoms to check, so no O2' found -> DNA logic
        return true;
    }
    
    bool has_o2_prime = false;
    for (const auto& atom : residue->atoms) {
        if (atom.name == " O2'" || atom.name == " O2*") {
            has_o2_prime = true;
            break;
        }
    }
    
    if (has_o2_prime)
        return false;
    else
        return true;
}

} // namespace util
} // namespace coot