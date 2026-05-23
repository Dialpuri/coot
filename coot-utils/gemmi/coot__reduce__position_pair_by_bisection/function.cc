#include "function.hh"
#include <clipper/core/coords.h>
#include <cmath>

namespace coot {

inline clipper::Coord_orth co_gemmi(const gemmi::Atom* at) {
    return clipper::Coord_orth(at->pos.x, at->pos.y, at->pos.z);
}

namespace reduce {

std::pair<clipper::Coord_orth, clipper::Coord_orth>
position_pair_by_bisection_gemmi(const gemmi::Atom* at_1,
                                  const gemmi::Atom* at_2,
                                  const gemmi::Atom* at_3,
                                  double bl,
                                  double alpha) {
    clipper::Coord_orth at_1_pos = co_gemmi(at_1);
    clipper::Coord_orth at_2_pos = co_gemmi(at_2);
    clipper::Coord_orth at_3_pos = co_gemmi(at_3);

    clipper::Coord_orth vec_1_to_2(at_2_pos - at_1_pos);
    clipper::Coord_orth vec_3_to_2(at_2_pos - at_3_pos);
    clipper::Coord_orth vec_1_to_2_uv(vec_1_to_2.unit());
    clipper::Coord_orth vec_3_to_2_uv(vec_3_to_2.unit());

    clipper::Coord_orth bisect_delta = 0.5 * (vec_1_to_2_uv + vec_3_to_2_uv);
    clipper::Coord_orth bisect_delta_uv(bisect_delta.unit());

    clipper::Coord_orth vec_1_to_3 = at_3_pos - at_1_pos;
    clipper::Coord_orth vec_1_to_3_uv(vec_1_to_3.unit());

    clipper::Coord_orth cpu(clipper::Coord_orth::cross(vec_1_to_3_uv, bisect_delta_uv));

    double scale_fac_bisector = bl * std::sin(0.5 * (M_PI - alpha));
    double scale_fac_cpu      = bl * std::cos(0.5 * (M_PI - alpha));

    clipper::Coord_orth Hp1 = at_2_pos + scale_fac_bisector * bisect_delta_uv - scale_fac_cpu * cpu;
    clipper::Coord_orth Hp2 = at_2_pos + scale_fac_bisector * bisect_delta_uv + scale_fac_cpu * cpu;

    return std::pair<clipper::Coord_orth, clipper::Coord_orth>(Hp1, Hp2);
}

} // namespace reduce
} // namespace coot