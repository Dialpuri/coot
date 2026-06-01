#pragma once

#include <vector>
#include <gemmi/model.hpp>

namespace coot {

struct restraints_container_t_gemmi {
    std::vector<gemmi::Atom> atom;
    std::vector<bool> fixed_during_refinement;
    size_t n_atoms = 0;

    void unset_fixed_during_refinement_udd_gemmi();
};

inline void restraints_container_t_gemmi::unset_fixed_during_refinement_udd_gemmi() {
    for (size_t i = 0; i < n_atoms; i++) {
        fixed_during_refinement[i] = false;
    }
}

} // namespace coot