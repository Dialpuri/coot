#pragma once
#include <gemmi/model.hpp>
#include <cmath>

namespace coot {

// Ported from MMDB version: distance(mmdb::Atom*, mmdb::Atom*)
// Calculates Euclidean distance between two atoms.
// Returns -1 if either atom is null.
inline double
distance_gemmi(const gemmi::Atom* at_1, const gemmi::Atom* at_2) {
    if (!at_1 || !at_2) {
        return -1.0;
    }
    double dx = at_1->pos.x - at_2->pos.x;
    double dy = at_1->pos.y - at_2->pos.y;
    double dz = at_1->pos.z - at_2->pos.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

}