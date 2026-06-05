#pragma once
#include <gemmi/model.hpp>
#include <string>

namespace coot {

// Port of coot::restraints_container_t::O3prime_and_P_are_close_ng
// Checks whether O3' atom in res_1 and P atom in res_2 are within d_crit distance.
// gemmi does not store per-atom TER flags (TER records cause chain breaks instead),
// so the TER check from the MMDB version is omitted.
inline bool O3prime_and_P_are_close_ng_gemmi(gemmi::Residue *res_1,
                                             gemmi::Residue *res_2,
                                             float d_crit) {
    bool status = false;

    // Find O3' atom in res_1
    const gemmi::Atom *at_1 = nullptr;
    if (res_1) {
        for (const auto &atom : res_1->atoms) {
            if (atom.name == " O3'") {
                at_1 = &atom;
                break;
            }
        }
    }

    // Find P atom in res_2
    const gemmi::Atom *at_2 = nullptr;
    if (res_2) {
        for (const auto &atom : res_2->atoms) {
            if (atom.name == " P  ") {
                at_2 = &atom;
                break;
            }
        }
    }

    if (at_1) {
        if (at_2) {
            float dx = at_1->pos.x - at_2->pos.x;
            float dy = at_1->pos.y - at_2->pos.y;
            float dz = at_1->pos.z - at_2->pos.z;
            float dd = dx * dx + dy * dy + dz * dz;
            if (dd < d_crit * d_crit)
                status = true;
        }
    }

    // NOTE: gemmi does not store per-atom TER flags (TER records cause chain
    // breaks instead). The original MMDB code checked for TER atoms in res_1
    // and would set status=false. In gemmi this is unnecessary because residues
    // at chain boundaries from TER records would be in separate chains.

    return status;
}

} // namespace coot