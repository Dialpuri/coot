#pragma once

#include <gemmi/model.hpp>

namespace Bond_lines_container {
namespace Bond_lines_container {

struct bonded_atom_pair_t {
    gemmi::Atom *atom_1;
    gemmi::Atom *atom_2;
    double residual;

    static bonded_atom_pair_t bonded_atom_pair_t_gemmi(gemmi::Atom *a1, gemmi::Atom *a2, double r) {
        return {a1, a2, r};
    }
};

} // namespace Bond_lines_container
} // namespace Bond_lines_container