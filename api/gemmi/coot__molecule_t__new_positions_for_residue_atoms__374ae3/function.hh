#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include "api/moved-atom.hh"
#include <algorithm>
#include <cctype>

namespace coot {

/// Trim leading and trailing whitespace from a string.
static std::string trim(const std::string& s) {
    auto start = s.begin();
    while (start != s.end() && std::isspace(static_cast<unsigned char>(*start)))
        ++start;
    auto end = s.end();
    do {
        --end;
    } while (std::distance(start, end) > 0 && std::isspace(static_cast<unsigned char>(*end)));
    return std::string(start, end + 1);
}

/// Gemmi-port of coot::molecule_t::new_positions_for_residue_atoms.
///
/// For each moved_atom_t in *moved_atoms*, finds a matching atom in *residue_p*
/// (by atom name and alternate-configuration) and updates its coordinates.
/// Returns the number of atoms whose positions were changed.
///
/// NOTE: make_backup has no gemmi port; the do_backup flag is accepted but
/// currently has no effect (the original caller was responsible for backups
/// anyway).
inline int new_positions_for_residue_atoms_gemmi(
    gemmi::Residue* residue_p,
    const std::vector<coot::api::moved_atom_t>& moved_atoms,
    bool /*do_backup*/) {

    int n_atoms_moved = 0;

    if (residue_p) {
        // do_backup: no gemmi port for make_backup — skip silently.
        // (The original caller handles backups at a higher level anyway.)

        for (const auto& mva : moved_atoms) {
            // Trim mva.atom_name (gemmi stores unpadded names like "CA"
            // but moved_atom_t may hold padded names like " CA ").
            std::string trimmed_mva_name = trim(mva.atom_name);

            for (auto& atom : residue_p->atoms) {
                // gemmi atom.name is already unpadded (e.g. "CA")
                if (atom.name != trimmed_mva_name) continue;

                // Normalize alt_conf: MMDB uses "" for "no alt loc" but gemmi
                // uses '\0' (null char). Treat empty string and null char as equivalent.
                std::string mva_alt = mva.alt_conf;
                char alt_char = atom.altloc;
                if (alt_char == '\0') {
                    if (!mva_alt.empty()) continue;  // gemmi has no alt, mva does — skip
                } else {
                    std::string alt_conf_str(1, alt_char);
                    if (alt_conf_str != mva_alt) continue;
                }

                atom.pos.x = static_cast<double>(mva.x);
                atom.pos.y = static_cast<double>(mva.y);
                atom.pos.z = static_cast<double>(mva.z);
                n_atoms_moved++;
            }
        }
    }

    return n_atoms_moved;
}

}  // namespace coot