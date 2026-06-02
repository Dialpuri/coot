#pragma once

#include <cmath>
#include <iosfwd>
#include <gemmi/math.hpp>

namespace coot {

// gemmi port of coot::coot_mat44::is_close_to_unit_matrix()
// Tests whether a gemmi::Transform is close to the identity (unit) matrix.
inline bool is_close_to_unit_matrix_transform_gemmi(const gemmi::Transform& t) {
    float sum_dist = 0.0f;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            float val;
            if (i < 3 && j < 3)
                val = static_cast<float>(t.mat.a[i][j]);
            else if (i < 3 && j == 3)
                val = static_cast<float>(t.vec.at(i));
            else if (i == 3 && j < 3)
                val = 0.0f;
            else  // i==3 && j==3  — last element of affine matrix is always 1
                val = 1.0f;

            if (i == j) {
                if (i == 3)
                    sum_dist += std::fabs(val);       // original: fabsf(m[3][3])
                else
                    sum_dist += std::fabs(val - 1.0f); // original: fabsf(m[i][i]-1)
            } else {
                sum_dist += std::fabs(val);           // original: fabsf(m[i][j])
            }
        }
    }
    return sum_dist < 0.001f;
}

// gemmi port of the coot::coot_mat44 stream insertion operator
// Writes a gemmi::Transform in the coot::mat_44 text format.
inline std::ostream& operator_stream_insert_gemmi(std::ostream& s, const gemmi::Transform& t) {
    s << "coot::mat_44\n";
    for (int i = 0; i < 3; i++) {
        s << "[ " << t.mat.a[i][0] << " "
                  << t.mat.a[i][1] << " "
                  << t.mat.a[i][2] << " "
                  << t.vec.at(i) << " ]\n";
    }
    s << "[ 0 0 0 1 ]\n";
    return s;
}

}  // namespace coot