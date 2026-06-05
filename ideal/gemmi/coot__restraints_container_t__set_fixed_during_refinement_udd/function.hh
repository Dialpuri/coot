#pragma once

#include <gemmi/model.hpp>
#include <set>
#include <vector>

namespace coot {

inline void set_fixed_during_refinement_udd_gemmi(
    std::vector<gemmi::Atom*>& atoms,
    int n_atoms,
    const std::set<int>& fixed_atom_indices,
    std::vector<bool>& fixed_during_refinement
) {
    fixed_during_refinement.assign(n_atoms, false);
    for (int i = 0; i < n_atoms; i++) {
        if (fixed_atom_indices.find(i) != fixed_atom_indices.end()) {
            fixed_during_refinement[i] = true;
        }
    }
}

} // namespace coot