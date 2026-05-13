#pragma once

#include <gemmi/model.hpp>
#include <clipper/clipper.h>

namespace coot { namespace util {

// Transform all atoms in a residue using a clipper coordinate transformation
// MMDB equivalent: coot::util::transform_atoms(mmdb::Residue *, const clipper::RTop_orth &)
inline void transform_atoms_gemmi(gemmi::Residue *res, const clipper::RTop_orth &rtop) {
    for (auto& atom : res->atoms) {
        clipper::Coord_orth co(atom.pos.x, atom.pos.y, atom.pos.z);
        clipper::Coord_orth trans_pos = co.transform(rtop);
        atom.pos.x = trans_pos.x();
        atom.pos.y = trans_pos.y();
        atom.pos.z = trans_pos.z();
    }
}

}} // namespace coot::util