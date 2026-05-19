#pragma once
#include <gemmi/model.hpp>

namespace coot {
namespace cho {

// Port of coot::cho::is_het_residue from MMDB to gemmi
// Returns true if the residue is a heterogen residue (het_flag == 'H')
inline bool is_het_residue_gemmi(const gemmi::Residue* res) {
    if (!res) return false;
    return res->het_flag == 'H';
}

} // namespace cho
} // namespace coot