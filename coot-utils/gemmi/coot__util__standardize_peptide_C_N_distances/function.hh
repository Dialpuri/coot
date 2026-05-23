#pragma once

#include <vector>
#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot {
namespace util {

// Port of coot::util::standardize_peptide_C_N_distances from MMDB to gemmi
// Standardizes C-N peptide bond distances to 1.320 Angstroms
inline void standardize_peptide_C_N_distances_gemmi(
    const std::vector<std::pair<gemmi::Atom*, gemmi::Atom*>>& C_N_pairs) {
    
    for (size_t i = 0; i < C_N_pairs.size(); i++) {
        gemmi::Atom* c_at = C_N_pairs[i].first;
        gemmi::Atom* n_at = C_N_pairs[i].second;
        if (c_at && n_at) {
            clipper::Coord_orth pt_1 = co_gemmi(c_at);
            clipper::Coord_orth pt_2 = co_gemmi(n_at);
            double d = clipper::Coord_orth::length(pt_1, pt_2);
            if (d > 1.32) {
                clipper::Coord_orth delta = pt_2 - pt_1;
                clipper::Coord_orth uv(delta.unit());
                double shift = 0.5 * (d - 1.32);
                c_at->pos.x += shift * uv.x();
                c_at->pos.y += shift * uv.y();
                c_at->pos.z += shift * uv.z();
                n_at->pos.x -= shift * uv.x();
                n_at->pos.y -= shift * uv.y();
                n_at->pos.z -= shift * uv.z();
            }
        }
    }
}

} // namespace util
} // namespace coot