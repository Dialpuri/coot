#pragma once

#include <gemmi/model.hpp>
#include <clipper/clipper.h>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot {

inline void transform_by_gemmi(gemmi::Structure& st, const clipper::RTop_orth &rtop) {
    for (gemmi::Model& model : st.models) {
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& residue : chain.residues) {
                for (gemmi::Atom& at : residue.atoms) {
                    // Skip termination markers (TMR atoms) - equivalent to !at->isTer()
                    if (at.name != "TMR") {
                        clipper::Coord_orth pos = coot::co_gemmi(&at);
                        clipper::Coord_orth p2 = pos.transform(rtop);
                        at.pos.x = p2.x();
                        at.pos.y = p2.y();
                        at.pos.z = p2.z();
                    }
                }
            }
        }
    }
}

} // namespace coot