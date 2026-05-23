#include "function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot {
namespace reduce {

clipper::Coord_orth position_by_bisection_gemmi(const gemmi::Atom* at_1,
                                                const gemmi::Atom* at_2,
                                                const gemmi::Atom* at_3,
                                                double bl) {
    clipper::Coord_orth at_1_pos = coot::co_gemmi(at_1);
    clipper::Coord_orth at_2_pos = coot::co_gemmi(at_2);
    clipper::Coord_orth at_3_pos = coot::co_gemmi(at_3);

    clipper::Coord_orth vec_1_to_2 = at_2_pos - at_1_pos;
    clipper::Coord_orth vec_3_to_2 = at_2_pos - at_3_pos;
    
    clipper::Coord_orth vec_1_to_2_uv = clipper::Coord_orth(vec_1_to_2.unit());
    clipper::Coord_orth vec_3_to_2_uv = clipper::Coord_orth(vec_3_to_2.unit());

    clipper::Coord_orth bisect_delta = 0.5 * (vec_1_to_2_uv + vec_3_to_2_uv);
    clipper::Coord_orth bisect_delta_uv = clipper::Coord_orth(bisect_delta.unit());

    clipper::Coord_orth Hp1 = at_2_pos + bl * bisect_delta_uv;

    return Hp1;
}

} // namespace reduce
} // namespace coot