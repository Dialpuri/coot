#pragma once

#include <gemmi/model.hpp>
#include <gemmi/math.hpp>

namespace coot {

inline double
angle_gemmi(const gemmi::Atom *at_1, const gemmi::Atom *at_2, const gemmi::Atom *at_3) {
    if (!at_1 || !at_2 || !at_3) {
        return -1;
    }

    gemmi::Vec3 pt_1 = at_1->pos;
    gemmi::Vec3 pt_2 = at_2->pos;
    gemmi::Vec3 pt_3 = at_3->pos;

    gemmi::Vec3 v1 = pt_1 - pt_2;
    gemmi::Vec3 v2 = pt_3 - pt_2;

    double angle_rad = v1.angle(v2);
    double angle_deg = angle_rad * 180.0 / M_PI;

    return angle_deg;
}

} // namespace coot