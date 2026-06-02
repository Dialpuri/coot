#pragma once

#include <gemmi/model.hpp>
#include <clipper/core/coords.h>

namespace coot {

inline void transform_by_gemmi(const clipper::RTop_orth &rtop, gemmi::Residue &residue_moving) {
    for (gemmi::Atom &atom : residue_moving.atoms) {
        clipper::Coord_orth p(atom.pos.x, atom.pos.y, atom.pos.z);
        clipper::Coord_orth p2 = p.transform(rtop);
        atom.pos.x = p2.x();
        atom.pos.y = p2.y();
        atom.pos.z = p2.z();
    }
}

} // namespace coot