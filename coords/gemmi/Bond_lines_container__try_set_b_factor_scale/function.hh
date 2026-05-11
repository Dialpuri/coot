#pragma once

#include <string>
#include <gemmi/pdb.hpp>

namespace coot {

class Bond_lines_container {
public:
    float b_factor_scale = 1.0f;
    
    void try_set_b_factor_scale_gemmi(const gemmi::Structure& st) {
        // Original MMDB version reads UDD data which gemmi doesn't support.
        // The function would read a scale factor from user-defined data,
        // but gemmi has no equivalent mechanism, so we keep the default.
        (void)st; // unused
    }
    
    static std::string b_factor_bonds_scale_handle_name;
};

} // namespace coot