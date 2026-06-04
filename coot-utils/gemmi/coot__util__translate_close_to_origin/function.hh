#pragma once
#include <gemmi/model.hpp>
#include <gemmi/unitcell.hpp>
#include <clipper/clipper.h>
#include <iostream>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_cell_symm/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__shift_to_origin/gemmi/function.hh"

namespace coot {
namespace util {

inline void translate_close_to_origin_gemmi(gemmi::Structure& st) {
    try {
        std::pair<clipper::Cell, clipper::Spacegroup> csp = get_cell_symm_gemmi(st);
        clipper::Coord_frac cf = shift_to_origin_gemmi(st);
        clipper::Coord_orth co = cf.coord_orth(csp.first);

        for (gemmi::Model& model : st.models) {
            for (gemmi::Chain& chain : model.chains) {
                for (gemmi::Residue& residue : chain.residues) {
                    for (gemmi::Atom& atom : residue.atoms) {
                        atom.pos.x += co.x();
                        atom.pos.y += co.y();
                        atom.pos.z += co.z();
                    }
                }
            }
        }
    } catch (const std::runtime_error& rte) {
        std::cout << rte.what() << std::endl;
    }
}

} // namespace util
} // namespace coot