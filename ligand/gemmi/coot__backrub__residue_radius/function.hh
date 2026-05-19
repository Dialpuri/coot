#pragma once

#include <string>
#include <vector>
#include <clipper/core/coords.h>
#include <cmath>

namespace coot { namespace backrub {

// Gemmi port of coot::backrub::residue_radius
// Takes a gemmi Residue instead of mmdb::Residue
inline float residue_radius_gemmi(const std::vector<clipper::Coord_orth>& residue_atoms,
                                  const clipper::Coord_orth& rc) {
    float longest_length = 0.0f;
    for (const auto& atom : residue_atoms) {
        float dx = atom.x() - rc.x();
        float dy = atom.y() - rc.y();
        float dz = atom.z() - rc.z();
        float this_length_sq = dx*dx + dy*dy + dz*dz;
        if (this_length_sq > longest_length) {
            longest_length = this_length_sq;
        }
    }
    return std::sqrt(longest_length);
}

}} // namespace coot::backrub