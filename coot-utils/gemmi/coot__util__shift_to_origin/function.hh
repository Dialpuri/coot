#pragma once
#include <clipper/clipper.h>
#include <gemmi/model.hpp>
#include <gemmi/unitcell.hpp>
#include <vector>
#include <algorithm>
#include <cmath>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_cell_symm/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__median_position/gemmi/function.hh"

namespace coot {
namespace util {

inline clipper::Coord_frac shift_to_origin_gemmi(const gemmi::Structure& st) {
    // Get cell and spacegroup (throws if no cell or symmetry)
    std::pair<clipper::Cell, clipper::Spacegroup> csp = get_cell_symm_gemmi(st);
    clipper::Cell cell = csp.first;

    // Collect all atom positions
    std::vector<gemmi::Position> positions;
    for (const auto& model : st.models) {
        for (const auto& chain : model.chains) {
            for (const auto& residue : chain.residues) {
                for (const auto& atom : residue.atoms) {
                    positions.push_back(atom.pos);
                }
            }
        }
    }

    // Compute median position
    gemmi::Position median_pos = median_position_gemmi(positions);

    // Convert to fractional coordinates using clipper::Cell
    clipper::Coord_orth orth(median_pos.x, median_pos.y, median_pos.z);
    clipper::Coord_frac mpf = orth.coord_frac(cell);

    // Round and negate
    clipper::Coord_frac rf(
        std::round(-mpf.u()),
        std::round(-mpf.v()),
        std::round(-mpf.w())
    );
    return rf;
}

} // namespace util
} // namespace coot