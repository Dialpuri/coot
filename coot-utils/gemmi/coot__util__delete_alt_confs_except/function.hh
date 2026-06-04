#pragma once

#include <gemmi/model.hpp>
#include <algorithm>
#include <string>

namespace coot {
namespace util {

// Delete atoms in *residue* whose altloc does not match *alt_conf*.
//
// MMDB convention: empty altLoc ("") means "no alternative conformation".
// gemmi convention: '\0' or ' ' means "no alternative conformation".
// Therefore an empty alt_conf string matches both '\0' and ' '.
inline void delete_alt_confs_except_gemmi(gemmi::Residue& residue,
                                          const std::string& alt_conf) {
  bool keep_empty = alt_conf.empty();
  char target = keep_empty ? '\0' : alt_conf[0];  // dummy when keep_empty
  auto& atoms = residue.atoms;
  atoms.erase(
      std::remove_if(atoms.begin(), atoms.end(),
                     [&target, keep_empty](const gemmi::Atom& a) {
                       if (keep_empty) {
                         // keep atoms with no alt-location
                         return !(a.altloc == '\0' || a.altloc == ' ');
                       }
                       return a.altloc != target;
                     }),
      atoms.end());
}

} // namespace util
} // namespace coot