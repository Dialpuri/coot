#pragma once

#include <gemmi/model.hpp>

namespace coot {
namespace util {

// Deep copy a residue (MMDB → gemmi port)
// Returns a newly allocated gemmi::Residue with all atoms copied
gemmi::Residue *deep_copy_this_residue_gemmi(const gemmi::Residue *residue) {
    if (!residue) {
        return nullptr;
    }
    
    gemmi::Residue *rres = new gemmi::Residue;
    rres->seqid.num.value = residue->seqid.num.value;
    rres->name = residue->name;
    rres->seqid.icode = residue->seqid.icode;
    
    // Copy all atoms
    for (const gemmi::Atom& atom : residue->atoms) {
        rres->atoms.push_back(atom);
    }
    
    return rres;
}

} // namespace util
} // namespace coot