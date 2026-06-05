#pragma once

#include <string>
#include <gemmi/model.hpp>
#include "coot/geometry/protein-geometry.hh"

namespace coot {

/// Gemmi-port of restraints_container_t::get_atom_index.
///
/// Iterates over atoms in a gemmi::Residue and returns the index of the first
/// atom whose name matches `atom_name_in` (exact string match, including PDB-style
/// padding).
///
/// Returns -2 when no match is found (matching original MMDB semantics).
///
/// NOTE: gemmi::atom.name is unpadded (e.g. "N"), while the MMDB API returned
/// 4-char padded names (e.g. " N  "). We pad the gemmi name using
/// coot::atom_id_mmdb_expand before comparison.
int get_atom_index_gemmi(const std::string &atom_name_in,
                         const gemmi::Residue &residue)
{
    int idx = -2;
    for (size_t i = 0; i < residue.atoms.size(); i++) {
        const auto &atom = residue.atoms[i];
        std::string padded_name = atom_id_mmdb_expand(atom.name);
        if (padded_name == atom_name_in) {
            idx = static_cast<int>(i);
            break;
        }
    }
    return idx;
}

} // namespace coot