#pragma once

#include <vector>

namespace coot {

// Returns the index of the hydrogen atom among index_1, index_2, index_3
// if exactly one is hydrogen; otherwise returns -1.
// is_hydrogen_flags[i] should be true if atom at index i is a hydrogen.
inline int get_chiral_hydrogen_index_gemmi(
    int index_1, int index_2, int index_3,
    const std::vector<bool>& is_hydrogen_flags
) {
    int idx_chiral = -1;
    int n_hydrogen = 0;
    if (is_hydrogen_flags[index_1]) { n_hydrogen++; idx_chiral = index_1; }
    if (is_hydrogen_flags[index_2]) { n_hydrogen++; idx_chiral = index_2; }
    if (is_hydrogen_flags[index_3]) { n_hydrogen++; idx_chiral = index_3; }

    if (n_hydrogen == 1) {
        return idx_chiral;
    } else {
        return -1;
    }
}

} // namespace coot