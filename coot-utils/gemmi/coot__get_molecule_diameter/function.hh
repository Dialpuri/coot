#pragma once
#include <vector>
#include <cmath>
#include <gemmi/model.hpp>
#include <coot/analysis/stats.hh>
#include <coot/utils/coot-utils.hh>

namespace coot {

float get_molecule_diameter_gemmi(const std::vector<gemmi::Atom>& atoms) {
    float f = -1;
    int n_max = static_cast<int>(atoms.size());

    stats::single s;
    for (unsigned int i = 0; i < 1000; i++) {
        float f1 = coot::util::random_f();
        float f2 = coot::util::random_f();
        int ff_1 = static_cast<float>(n_max) * f1;
        int ff_2 = static_cast<float>(n_max) * f2;

        if (f1 < 1.0) {
            if (f2 < 1.0) {
                int idx_1 = static_cast<int>(ff_1);
                int idx_2 = static_cast<int>(ff_2);
                if (idx_1 != idx_2) {
                    const gemmi::Atom& at_1 = atoms[idx_1];
                    const gemmi::Atom& at_2 = atoms[idx_2];
                    float dx = at_2.pos.x - at_1.pos.x;
                    float dy = at_2.pos.y - at_1.pos.y;
                    float dz = at_2.pos.z - at_1.pos.z;
                    float dd = dx * dx + dy * dy + dz * dz;
                    float d = std::sqrt(dd);
                    s.add(d);
                }
            }
        }
    }

    if (s.size() > 10) {
        f = s.get_ith_highest(10);
    }

    return f;
}

}