#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include "clipper/core/coords.h"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_residue_alt_confs/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__position_by_tetrahedron/gemmi/function.hh"

namespace coot {
namespace reduce {

/// Find atom by name in a residue using relaxed matching (substring search)
/// plus exact altloc match.
static inline const gemmi::Atom* find_atom_ptr(
    const gemmi::Residue& residue, const std::string& name_str,
    char wanted_alt) {
    for (const auto& atom : residue.atoms) {
        if (atom.altloc != wanted_alt)
            continue;
        bool matched_name = (atom.name.size() <= 4)
            ? (atom.name == name_str
               || atom.name.substr(0, name_str.size()) == name_str
               || atom.name.find(name_str) != std::string::npos)
            : false;
        if (matched_name)
            return &atom;
    }
    return nullptr;
}

/// Add a hydrogen atom to a residue. Always appends a new atom,
/// replicating the MMDB behaviour where name padding causes the
/// existing-atom lookup to fail.
static inline gemmi::Atom* add_h_atom(
    const std::string& atom_name, const clipper::Coord_orth& pos,
    float bf, char altloc, gemmi::Residue& residue) {

    gemmi::Atom new_H;
    new_H.name = atom_name;
    new_H.element = gemmi::Element("H");
    new_H.pos = gemmi::Position(pos.x(), pos.y(), pos.z());
    new_H.occ = 1.0f;
    new_H.b_iso = bf;
    new_H.altloc = altloc;

    residue.atoms.push_back(new_H);
    return &residue.atoms.back();
}

/// 5-parameter overload: position hydrogen tetrahedrally.
inline void
add_tetrahedral_hydrogen_gemmi(const std::string& H_at_name,
                               const std::string& at_central_name,
                               const std::string& neighb_at_name_1,
                               const std::string& neighb_at_name_2,
                               const std::string& neighb_at_name_3,
                               double bond_length,
                               gemmi::Residue& residue) {

   std::vector<std::string> alt_confs = coot::util::get_residue_alt_confs_gemmi(residue);

   for (size_t i = 0; i < alt_confs.size(); i++) {
      char wanted_alt = alt_confs[i].empty() ? '\0' : alt_confs[i][0];

      const gemmi::Atom* at_central = find_atom_ptr(residue, at_central_name, wanted_alt);
      const gemmi::Atom* at_n_1 = find_atom_ptr(residue, neighb_at_name_1, wanted_alt);
      const gemmi::Atom* at_n_2 = find_atom_ptr(residue, neighb_at_name_2, wanted_alt);
      const gemmi::Atom* at_n_3 = find_atom_ptr(residue, neighb_at_name_3, wanted_alt);

      if (at_central && at_n_1 && at_n_2 && at_n_3) {
         clipper::Coord_orth H_pos = coot::reduce::position_by_tetrahedron_gemmi(
             at_central, at_n_1, at_n_2, at_n_3, bond_length);
         float bf = static_cast<float>(at_central->b_iso);
         add_h_atom(H_at_name, H_pos, bf, wanted_alt, residue);
      }
   }
}

/// Vector overload: only adds when exactly 3 neighbours are supplied.
inline void
add_tetrahedral_hydrogen_gemmi(const std::string& H_at_name,
                               const std::string& first_neighb,
                               const std::vector<std::string>& second_neighb_vec,
                               double bl,
                               gemmi::Residue& residue) {

   if (second_neighb_vec.size() == 3)
      add_tetrahedral_hydrogen_gemmi(H_at_name, first_neighb,
                                     second_neighb_vec[0],
                                     second_neighb_vec[1],
                                     second_neighb_vec[2],
                                     bl, residue);
   else
      std::cout << "WARNING:: atom " << first_neighb << " had " << second_neighb_vec.size()
                << " neighbours  (not 3)" << std::endl;
}

} // namespace reduce
} // namespace coot