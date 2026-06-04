#pragma once
#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot {
namespace util {

// Port of coot::util::transform_chain from MMDB to gemmi
// Applies an orthonormal rigid transformation to every atom in a chain
inline void transform_chain_gemmi(gemmi::Chain &chain, const clipper::RTop_orth &rtop) {
    for (gemmi::Residue &residue : chain.residues) {
        for (gemmi::Atom &at : residue.atoms) {
            clipper::Coord_orth pt(co_gemmi(&at));
            clipper::Coord_orth new_pt(rtop * pt);
            at.pos.x = new_pt.x();
            at.pos.y = new_pt.y();
            at.pos.z = new_pt.z();
        }
    }
}

} // namespace util
} // namespace coot