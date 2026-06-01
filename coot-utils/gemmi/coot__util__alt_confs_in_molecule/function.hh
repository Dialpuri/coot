#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <set>

namespace coot {
namespace util {

/**
 * Return the set of unique alternate-loc indicators found in the first
 * model of a gemmi::Structure.
 *
 * Equivalent to coot::util::alt_confs_in_molecule (MMDB version) but uses
 * gemmi types throughout.
 *
 * Mapping of no-alt-loc: gemmi stores ' ' for "no alternate location" while
 * MMDB stored '\0'. This function normalises ' ' → "" so that callers
 * expecting the MMDB convention get the same result.
 */
std::vector<std::string>
alt_confs_in_molecule_gemmi(gemmi::Structure const& st) {

   std::set<std::string> s;

   // The original MMDB version accessed model 1 (index 0 in gemmi).
   if (!st.models.empty()) {
      gemmi::Model const& model = st.models[0];
      for (gemmi::Chain const& chain : model.chains) {
         for (gemmi::Residue const& residue : chain.residues) {
            for (gemmi::Atom const& atom : residue.atoms) {
               // In gemmi, TER records are not atoms, so every entry in
               // residue.atoms is a real atom (no isTer() check needed).
               char loc = atom.altloc;
               // Normalise: gemmi may use ' ' or '\0' for "no alt loc" (MMDB
               // used '\0'). Convert both → "" so callers get the same
               // empty-string semantics as the MMDB version.
               std::string alt_conf;
               if (loc != ' ' && loc != '\0') {
                  alt_conf = std::string(1, loc);
               }
               s.insert(alt_conf);
            }
         }
      }
   }

   // Convert sorted set to vector (same as original).
   std::vector<std::string> v;
   for (std::set<std::string>::const_iterator it = s.begin();
        it != s.end(); ++it) {
      v.push_back(*it);
   }
   return v;
}

} // namespace util
} // namespace coot