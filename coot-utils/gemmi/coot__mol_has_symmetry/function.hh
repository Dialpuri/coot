#pragma once
#include <gemmi/pdb.hpp>
#include <string>

namespace coot {

// Check if a structure has crystal symmetry (space group information)
// Ported from mmdb::Manager::GetTMatrix approach
inline bool mol_has_symmetry_gemmi(const gemmi::Structure& st) {
    // In gemmi, spacegroup_hm is non-empty if the structure has symmetry
    return !st.spacegroup_hm.empty();
}

} // namespace coot