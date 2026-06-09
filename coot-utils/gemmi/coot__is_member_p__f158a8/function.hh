#pragma once

#include <gemmi/model.hpp>
#include <vector>

namespace coot {

/// Check whether a gemmi::Residue* is present in a vector of gemmi::Residue*
/// (gemmi port of coot::is_member_p)
inline bool is_member_p_gemmi(const std::vector<gemmi::Residue*> &v, gemmi::Residue* a) {
    for (unsigned int i = 0; i < v.size(); i++) {
        if (v[i] == a) {
            return true;
        }
    }
    return false;
}

} // namespace coot
