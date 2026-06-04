#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_residue_alt_confs/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__position_by_tetrahedron/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__add_hydrogen_atom/gemmi/function.hh"

namespace coot::reduce {

// Find an atom in a gemmi::Residue by its unpadded name and altloc.
// gemmi stores atom names WITHOUT padding (e.g. "CA" not " CA ").
// altconf: empty string or " " means match any atom with no altloc ('\0' or ' ').
inline const gemmi::Atom * find_atom_in_residue(const gemmi::Residue &residue,
                                                 const std::string &atom_name,
                                                 const std::string &altconf) {
  for (const auto &atom : residue.atoms) {
    if (atom.name == atom_name) {
      if (altconf.empty() || altconf == " ") {
        if (atom.altloc == '\0' || atom.altloc == ' ')
          return &atom;
      } else {
        if (atom.altloc == static_cast<char>(altconf[0]))
          return &atom;
      }
    }
  }
  return nullptr;
}

inline void add_main_chain_HA_gemmi(gemmi::Residue &residue) {
  double bl = 0.97;

  std::vector<std::string> alt_confs = coot::util::get_residue_alt_confs_gemmi(residue);

  // If residue has no alternative conformations (no altloc set on any atom),
  // get_residue_alt_confs_gemmi returns an empty vector. In that case we use
  // a single empty-string entry so the loop runs once.
  if (alt_confs.empty()) {
    alt_confs.push_back("");
  }

  for (unsigned int i = 0; i < alt_confs.size(); i++) {
    // gemmi atom names are unpadded: "CA", "C", "N", "CB"
    const gemmi::Atom *at_ca  = find_atom_in_residue(residue, "CA",  alt_confs[i]);
    const gemmi::Atom *at_n1  = find_atom_in_residue(residue, "C",   alt_confs[i]);
    const gemmi::Atom *at_n2  = find_atom_in_residue(residue, "N",   alt_confs[i]);
    const gemmi::Atom *at_n3  = find_atom_in_residue(residue, "CB",  alt_confs[i]);

    if (at_ca && at_n1 && at_n2 && at_n3) {
      clipper::Coord_orth pos = coot::reduce::position_by_tetrahedron_gemmi(
          at_ca, at_n1, at_n2, at_n3, bl);
      float bf = static_cast<float>(at_ca->b_iso);
      coot::reduce::add_hydrogen_atom_gemmi(" HA ", pos, bf, alt_confs[i], residue);
    }
  }
}

} // namespace coot::reduce