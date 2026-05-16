#pragma once

#include <vector>
#include <gemmi/math.hpp>

namespace coot {

inline
std::vector<gemmi::Vec3>
flip_internal_gemmi(const std::vector<gemmi::Vec3> &ca_in,
                    const std::vector<gemmi::Vec3> &atoms) {

    std::vector<gemmi::Vec3> atoms_orth(atoms.size());
    std::vector<gemmi::Vec3> cas = ca_in;

    gemmi::Vec3 trans = cas[0];

    cas[0] -= trans;
    cas[1] -= trans;

    for (size_t i = 0; i < atoms.size(); ++i) {
        atoms_orth[i] = atoms[i];
        atoms_orth[i] -= trans;
    }

    gemmi::Vec3 ca_vec_unit = cas[1];
    double len = cas[1].length();
    if (len > 0) {
        ca_vec_unit /= len;
    }

    double l = ca_vec_unit.x;
    double m = ca_vec_unit.y;
    double n = ca_vec_unit.z;

    double ll = l * l;
    double mm = m * m;
    double nn = n * n;

    // Rotation matrix for 180° around the CA-C axis (Householder-like)
    gemmi::Mat33 r(
        ll - (mm + nn),   2.0 * l * m,            2.0 * n * l,
        2.0 * l * m,      mm - (ll + nn),         2.0 * m * n,
        2.0 * n * l,      2.0 * m * n,            nn - (ll + mm)
    );

    for (size_t i = 0; i < atoms_orth.size(); ++i) {
        atoms_orth[i] = r.multiply(atoms_orth[i]);
        atoms_orth[i] += trans;
    }

    return atoms_orth;
}

}  // namespace coot