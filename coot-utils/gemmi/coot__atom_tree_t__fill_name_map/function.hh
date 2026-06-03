#pragma once

#include <gemmi/model.hpp>
#include <map>
#include <string>
#include "coot-utils/map-index.hh"

namespace coot {

// Gemmi port of coot::atom_tree_t::fill_name_map
// Iterates over atoms in a gemmi::Residue, filters by altconf,
// and builds a name_to_index map (atom name -> map_index_t).
std::map<std::string, map_index_t> fill_name_map_gemmi(
    const gemmi::Residue& res,
    const std::string& altconf)
{
    std::map<std::string, map_index_t> name_to_index;

    int iat = 0;
    for (const gemmi::Atom& atom : res.atoms) {
        // gemmi stores trimmed names; MMDB uses 4-char padded names.
        // padded_name() left-pads; right-pad to 4 to match MMDB format.
        std::string atom_name = atom.padded_name();
        while (atom_name.size() < 4) atom_name += ' ';

        char atom_altl = atom.altloc;

        // In gemmi, no altloc is '\0' or ' '.
        // Original MMDB logic: atom_altl == "" || atom_altl == altconf
        bool no_alt = (atom_altl == '\0' || atom_altl == ' ');
        if (no_alt || (!altconf.empty() && atom_altl == altconf[0])) {
            name_to_index[atom_name] = map_index_t(iat);
        }
        iat++;
    }

    return name_to_index;
}

} // namespace coot