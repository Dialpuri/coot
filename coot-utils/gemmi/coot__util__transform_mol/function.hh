#pragma once

#include <gemmi/model.hpp>
#include <clipper/clipper.h>

namespace coot { namespace util {

// MMDB equivalent: coot::util::transform_mol(mmdb::Manager *, const clipper::RTop_orth &)
inline void transform_mol_gemmi(gemmi::Structure *st, const clipper::RTop_orth &rtop) {
    for (auto& model : st->models) {
        for (auto& chain : model.chains) {
            for (auto& residue : chain.residues) {
                for (auto& atom : residue.atoms) {
                    clipper::Coord_orth co(atom.pos.x, atom.pos.y, atom.pos.z);
                    clipper::Coord_orth trans_pos = co.transform(rtop);
                    atom.pos.x = trans_pos.x();
                    atom.pos.y = trans_pos.y();
                    atom.pos.z = trans_pos.z();
                }
            }
        }
    }
}

}} // namespace coot::util