#pragma once
#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include <cstring>

namespace coot { namespace util {

// Deep copy a residue (returning a new Residue* that caller must delete)
// Returns nullptr if input is nullptr
inline gemmi::Residue*
deep_copy_this_residue_gemmi(const gemmi::Residue* residue) {
    if (!residue) {
        return nullptr;
    }

    gemmi::Residue* rres = new gemmi::Residue;
    rres->seqid.num.value = residue->seqid.num.value;
    rres->name = residue->name;
    rres->seqid.icode = residue->seqid.icode;

    for (const gemmi::Atom& atom : residue->atoms) {
        // In gemmi, there's no is_ter() - we check element instead
        // TER atoms have element "X" in gemmi's representation
        if (atom.element != gemmi::Element("X")) {
            gemmi::Atom atom_p = atom;
            rres->atoms.push_back(atom_p);
        }
    }

    return rres;
}

}} // namespace coot::util