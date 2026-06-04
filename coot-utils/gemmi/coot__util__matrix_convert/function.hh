#pragma once

#include <gemmi/math.hpp>
#include "clipper/core/coords.h"

namespace coot {
namespace util {

inline clipper::RTop_orth matrix_convert_gemmi(gemmi::Transform transform) {
    gemmi::Mat33 mat = transform.mat;
    gemmi::Vec3 trans = transform.vec;

    clipper::Mat33<double> clipper_mat(
        mat[0][0], mat[0][1], mat[0][2],
        mat[1][0], mat[1][1], mat[1][2],
        mat[2][0], mat[2][1], mat[2][2]);

    clipper::Coord_orth cco(trans.x, trans.y, trans.z);
    return clipper::RTop_orth(clipper_mat, cco);
}

} // namespace util
} // namespace coot