#pragma once

#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include "clipper/core/coords.h"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_residue_alt_confs/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__position_by_tetrahedron/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__add_hydrogen_atom/gemmi/function.hh"

// Helper to trim whitespace from both ends of a string
inline std::string trim_str(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t\n\r");
    return str.substr(first, last - first + 1);
}

namespace coot {
namespace reduce {

/// Add a tetrahedral hydrogen to a residue.
/// Given a central atom name and three neighbor atom names, positions
/// the hydrogen tetrahedrally relative to the central atom.
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
      const gemmi::Atom* at_central = nullptr;
      const gemmi::Atom* at_n_1 = nullptr;
      const gemmi::Atom* at_n_2 = nullptr;
      const gemmi::Atom* at_n_3 = nullptr;

      // alt_confs[i] is "" for default (no alt loc), otherwise "A", "B", etc.
      auto is_match_altloc = [&](char atom_alt) -> bool {
         if (alt_confs[i].empty()) {
            // default altloc: match ' ' or '\0'
            return (atom_alt == ' ' || atom_alt == '\0');
         }
         return atom_alt == alt_confs[i][0];
      };
      for (const auto& atom : residue.atoms) {
         if (!is_match_altloc(atom.altloc)) continue;
         std::string trimmed = trim_str(atom.name);
         if (trimmed == trim_str(at_central_name))  at_central = &atom;
         if (trimmed == trim_str(neighb_at_name_1)) at_n_1 = &atom;
         if (trimmed == trim_str(neighb_at_name_2)) at_n_2 = &atom;
         if (trimmed == trim_str(neighb_at_name_3)) at_n_3 = &atom;
      }

      if (at_central && at_n_1 && at_n_2 && at_n_3) {
         clipper::Coord_orth H_pos = coot::reduce::position_by_tetrahedron_gemmi(
             at_central, at_n_1, at_n_2, at_n_3, bond_length);
         float bf = static_cast<float>(at_central->b_iso);
         coot::reduce::add_hydrogen_atom_gemmi(H_at_name, H_pos, bf, alt_confs[i], residue);
      }
   }
}

} // namespace reduce
} // namespace coot