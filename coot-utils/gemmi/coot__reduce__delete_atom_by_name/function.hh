#pragma once

#include <string>
#include <gemmi/model.hpp>

namespace coot {
namespace reduce {

// Delete all hydrogen atoms (element "H" or "D") from the given residue.
// Repeatedly scans and deletes until no more hydrogens are found.
inline void delete_atom_by_name_gemmi(const std::string & /* at_name */, gemmi::Residue &residue) {
  bool an_atom_was_deleted = true;
  while (an_atom_was_deleted) {
    an_atom_was_deleted = false;
    for (auto it = residue.atoms.begin(); it != residue.atoms.end(); ) {
      std::string ele = it->element.name();
      if (ele == "H" || ele == "D") {
        it = residue.atoms.erase(it);
        an_atom_was_deleted = true;
      } else {
        ++it;
      }
    }
  }
}

} // namespace reduce
} // namespace coot