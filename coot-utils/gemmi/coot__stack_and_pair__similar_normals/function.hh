#pragma once
#include <map>
#include <gemmi/model.hpp>
#include <clipper/core/coords.h>

namespace coot {
namespace stack_and_pair {

// Port of similar_normals from mmdb to gemmi
// Uses gemmi::Residue* instead of mmdb::Residue*
// Uses clipper::Coord_orth (same as original, since gemmi uses clipper for some types)
// The function is part of coot::stack_and_pair class, so we need the angle_crit member
inline bool similar_normals_gemmi(gemmi::Residue *res_1, gemmi::Residue *res_2,
                                  const std::map<gemmi::Residue *, clipper::Coord_orth> &normal_map,
                                  double angle_crit) {
    // normals are presumed to be normalized.
    bool status = false;
    std::map<gemmi::Residue *, clipper::Coord_orth>::const_iterator it_1;
    std::map<gemmi::Residue *, clipper::Coord_orth>::const_iterator it_2;

    it_1 = normal_map.find(res_1);
    it_2 = normal_map.find(res_2);

    if (it_1 != normal_map.end()) {
        if (it_2 != normal_map.end()) {
            const clipper::Coord_orth &n1 = it_1->second;
            const clipper::Coord_orth &n2 = it_2->second;
            double dp = clipper::Coord_orth::dot(n1, n2);
            double cos_angle_crit = cos(angle_crit);
            if ((dp > cos_angle_crit) || (dp < -cos_angle_crit)) {
                status = true;
            }
        }
    }
    return status;
}

} // namespace stack_and_pair
} // namespace coot