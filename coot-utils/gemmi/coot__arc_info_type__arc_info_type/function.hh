#pragma once
#include <gemmi/model.hpp>
#include <gemmi/math.hpp>
#include <stdexcept>
#include <cmath>

namespace coot {

struct arc_info_type {
    gemmi::Vec3 normal;
    gemmi::Vec3 start_point;
    gemmi::Vec3 start_dir;
    double delta;
    gemmi::Mat33 orientation_matrix;

    arc_info_type(const gemmi::Atom* at_1, const gemmi::Atom* at_2, const gemmi::Atom* at_3) {
        if (!at_1) throw std::runtime_error("null at_1");
        if (!at_2) throw std::runtime_error("null at_2");
        if (!at_3) throw std::runtime_error("null at_3");

        gemmi::Vec3 p1(at_1->pos.x, at_1->pos.y, at_1->pos.z);
        gemmi::Vec3 p2(at_2->pos.x, at_2->pos.y, at_2->pos.z);
        gemmi::Vec3 p3(at_3->pos.x, at_3->pos.y, at_3->pos.z);

        gemmi::Vec3 v1((p3 - p2).normalized());
        gemmi::Vec3 v2((p1 - p2).normalized());

        gemmi::Vec3 v3 = (v1.cross(v2)).normalized();
        gemmi::Vec3 v4 = (v3.cross(v1)).normalized();

        normal = v3;
        start_point = p2;
        start_dir = v2.normalized();

        gemmi::Vec3 d1 = (p1 - p2).normalized();
        gemmi::Vec3 d2 = (p3 - p2).normalized();
        double dot = d1.dot(d2);
        double angle_rad = std::acos(std::clamp(dot, -1.0, 1.0));
        delta = angle_rad * 180.0 / M_PI;

        orientation_matrix = gemmi::Mat33{
            v1.x, v1.y, v1.z,
            v4.x, v4.y, v4.z,
            v3.x, v3.y, v3.z
        };
    }
};

} // namespace coot