#pragma once

#include <string>
#include <gemmi/model.hpp>

namespace coot {

/// Gemmi-port of restraints_container_t::get_atom_index_for_restraint_using_alt_conf.
///
/// Iterates over atoms in a gemmi::Residue and returns the index of the first
/// atom whose name matches `atom_name` (exact string match, including PDB-style
/// padding) and whose alt_conf matches `alt_conf`.
///
/// alt_conf rules (mirroring MMDB semantics):
///   - An atom with no alt_conf (gemmi: altloc == ' ') always matches.
///   - An atom with alt_conf matches only when alt_conf[0] equals the atom's altloc.
///
/// Returns -1 when no match is found.
///
/// NOTE: MMDB's GetUDData (user-defined data) has no gemmi equivalent. This
/// function returns the local atom index within the residue's atom vector
/// instead.
int get_atom_index_for_restraint_using_alt_conf_gemmi(
    const std::string &atom_name,
    const std::string &alt_conf,
    const gemmi::Residue &residue)
{
    int idx = -1;
    for (size_t i = 0; i < residue.atoms.size(); i++) {
        const auto &atom = residue.atoms[i];
        if (atom.name == atom_name) {
            // Match alt_conf:
            //   - atom with no alt_conf (altloc == ' ') → matches anything
            //   - atom with alt_conf → matches only if alt_conf[0] == atom.altloc
            bool alt_matches = (atom.altloc == ' ') ||
                               (!alt_conf.empty() && atom.altloc == alt_conf[0]);
            if (alt_matches) {
                idx = static_cast<int>(i);
            }
        }
    }
    return idx;
}

} // namespace coot