#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <utility>
#include "coot-utils/atom-tree.hh"
#include "geometry/protein-geometry.hh"

namespace coot {

inline std::pair<bool, atom_index_quad>
get_atom_index_quad_gemmi(
    const dict_torsion_restraint_t &tr,
    gemmi::Residue &res,
    const std::string &altconf) {
  atom_index_quad quad(-1, -1, -1, -1);

  for (int iat = 0; iat < static_cast<int>(res.atoms.size()); iat++) {
    const auto &atom = res.atoms[iat];

    // altloc match: match if no altconf specified (accept any),
    // or atom has no altloc ('\0' or ' '), or it matches requested altconf
    bool this_alt_match = altconf.empty();
    if (!altconf.empty()) {
      char atom_altloc = atom.altloc;
      this_alt_match = (atom_altloc == '\0' || atom_altloc == ' ' ||
                        atom_altloc == static_cast<char>(altconf[0]));
    }

    if (!this_alt_match) continue;

    // gemmi atom.name is trimmed; _4c() returns 4-char padded.
    // Use atom_id_mmdb_expand to bring gemmi names to the same format.
    std::string padded_name = atom_id_mmdb_expand(atom.name);
    if (padded_name == tr.atom_id_1_4c()) quad.index1 = iat;
    if (padded_name == tr.atom_id_2_4c()) quad.index2 = iat;
    if (padded_name == tr.atom_id_3_4c()) quad.index3 = iat;
    if (padded_name == tr.atom_id_4_4c()) quad.index4 = iat;
  }

  bool success = (quad.index1 != -1 && quad.index2 != -1 &&
                  quad.index3 != -1 && quad.index4 != -1);
  return std::make_pair(success, quad);
}

} // namespace coot