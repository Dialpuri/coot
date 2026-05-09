#pragma once
#include <gemmi/model.hpp>

namespace coot {
namespace molecule_t {

// Port of mmdb::Residue* version to gemmi
// Returns true if the residue is a heterogen residue (het_flag == 'H')
inline bool is_het_residue_gemmi(const gemmi::Residue* res) {
    if (!res) return false;
    return res->het_flag == 'H';
}

// Alternative that takes a Residue reference (more idiomatic gemmi)
inline bool is_het_residue_gemmi(const gemmi::Residue& res) {
    return res.het_flag == 'H';
}

} // namespace molecule_t
} // namespace coot