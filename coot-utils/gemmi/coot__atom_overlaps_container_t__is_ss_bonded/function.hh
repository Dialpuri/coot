#pragma once
#include <gemmi/model.hpp>

namespace coot {

// gemmi port of atom_overlaps_container_t::is_ss_bonded
// Returns true if the residue is CYS (potential disulfide-bonded residue),
// false otherwise (including nullptr).
bool is_ss_bonded_gemmi(const gemmi::Residue* residue_p) {
    bool status = false;
    if (residue_p) {
        std::string res_name = residue_p->name;
        if (res_name == "CYS") {
            status = true;
        }
    }
    return status;
}

} // namespace coot