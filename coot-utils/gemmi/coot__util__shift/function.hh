#pragma once

#include <gemmi/model.hpp>
#include <clipper/core/coords.h>

namespace coot {
namespace util {

// Shift all atom coordinates in all models of `st` by `shift`.
inline void
shift_gemmi(gemmi::Structure &st, clipper::Coord_orth shift) {
    for (auto &model : st.models) {
        for (auto &chain : model.chains) {
            for (auto &residue : chain.residues) {
                for (auto &atom : residue.atoms) {
                    atom.pos.x += shift.x();
                    atom.pos.y += shift.y();
                    atom.pos.z += shift.z();
                }
            }
        }
    }
}

} // namespace util
} // namespace coot