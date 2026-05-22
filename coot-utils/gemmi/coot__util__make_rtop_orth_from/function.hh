#pragma once
#include <clipper/core/coords.h>
#include <clipper/core/clipper_types.h>

namespace coot { namespace util {

// Ported from MMDB version: converts 4x4 matrix to clipper::RTop_orth
// Uses plain C++ array instead of mmdb::mat44
clipper::RTop_orth
make_rtop_orth_from_gemmi(double mat[4][4]) {

   clipper::Mat33<double> clipper_mat(mat[0][0], mat[0][1], mat[0][2],
                                      mat[1][0], mat[1][1], mat[1][2],
                                      mat[2][0], mat[2][1], mat[2][2]);
   clipper::Coord_orth  cco(mat[0][3], mat[1][3], mat[2][3]);
   clipper::RTop_orth rtop(clipper_mat, cco);

   return rtop;
}

}} // namespace coot::util