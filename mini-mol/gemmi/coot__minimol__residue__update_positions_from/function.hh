#pragma once

#include <gemmi/model.hpp>
#include "mini-mol/mini-mol.hh"
#include "clipper/core/coords.h"

namespace coot { namespace minimol {

// Free-function gemmi port of the original member function
// residue::update_positions_from(mmdb::Residue *residue_p)
inline void update_positions_from_gemmi(residue &minires, const gemmi::Residue *residue_p) {
    int n_atoms = static_cast<int>(residue_p->atoms.size());
    // presuming the same atom order (don't check atom names) - funny
    // things will happen if this is not true.
    if (static_cast<int>(minires.atoms.size()) == n_atoms) {
        for (int iat = 0; iat < n_atoms; iat++) {
            const gemmi::Atom &at = residue_p->atoms[iat];
            clipper::Coord_orth p(at.pos.x, at.pos.y, at.pos.z);
            minires.atoms[iat].pos = p;
        }
    }
}

}} // namespace coot::minimol