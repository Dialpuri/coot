#pragma once

#include <gemmi/model.hpp>
#include <vector>

namespace coot {

struct torsion_atom_quad_gemmi {
    gemmi::Atom *atom_1;
    gemmi::Atom *atom_2;
    gemmi::Atom *atom_3;
    gemmi::Atom *atom_4;

    torsion_atom_quad_gemmi()
        : atom_1(nullptr), atom_2(nullptr), atom_3(nullptr), atom_4(nullptr) {}

    torsion_atom_quad_gemmi(gemmi::Atom *a1, gemmi::Atom *a2, gemmi::Atom *a3, gemmi::Atom *a4,
                            double angle_in, double angle_esd_in, int period_in)
        : atom_1(a1), atom_2(a2), atom_3(a3), atom_4(a4),
          angle(angle_in), angle_esd(angle_esd_in), period(period_in) {}

    double angle;
    double angle_esd;
    int period;
};

// Check whether both atoms appear together in any torsion quad.
// This is the gemmi-port of coot::both_in_a_torsion_p.
inline bool both_in_a_torsion_p_gemmi(gemmi::Atom *at_1,
                                      gemmi::Atom *at_2,
                                      const std::vector<torsion_atom_quad_gemmi> &quads) {
    bool in_a_tors = false;
    for (unsigned int i = 0; i < quads.size(); i++) {
        bool found_at_1 = false;
        bool found_at_2 = false;
        const torsion_atom_quad_gemmi &q = quads[i];
        if (q.atom_1 == at_1) found_at_1 = true;
        if (q.atom_2 == at_1) found_at_1 = true;
        if (q.atom_3 == at_1) found_at_1 = true;
        if (q.atom_4 == at_1) found_at_1 = true;
        if (q.atom_1 == at_2) found_at_2 = true;
        if (q.atom_2 == at_2) found_at_2 = true;
        if (q.atom_3 == at_2) found_at_2 = true;
        if (q.atom_4 == at_2) found_at_2 = true;
        if (found_at_1 && found_at_2) {
            in_a_tors = true;
            break;
        }
    }
    return in_a_tors;
}

} // namespace coot