#pragma once
#include <utility>
#include <gemmi/model.hpp>
#include <clipper/core/coords.h>

namespace coot {

// Calculate the closest approach distance between two residues.
// Returns (success_flag, distance_in_angstroms).
// If either residue has no atoms, returns (false, 0.0).
std::pair<bool, float> closest_approach_gemmi(const gemmi::Residue& r1,
                                               const gemmi::Residue& r2) {
    if (r1.atoms.empty() || r2.atoms.empty())
        return std::make_pair(false, 0.0f);
    
    float dist_sqrd_best = 9999999.9f;
    bool good_d = false;
    
    for (const auto& atom1 : r1.atoms) {
        for (const auto& atom2 : r2.atoms) {
            clipper::Coord_orth pt1(atom1.pos.x, atom1.pos.y, atom1.pos.z);
            clipper::Coord_orth pt2(atom2.pos.x, atom2.pos.y, atom2.pos.z);
            double d_sqrd = (pt2 - pt1).lengthsq();
            
            if (d_sqrd < dist_sqrd_best) {
                dist_sqrd_best = static_cast<float>(d_sqrd);
                good_d = true;
            }
        }
    }
    
    return std::make_pair(good_d, std::sqrt(dist_sqrd_best));
}

}