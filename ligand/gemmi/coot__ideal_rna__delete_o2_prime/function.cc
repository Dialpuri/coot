#include "function.hh"
#include <algorithm>

namespace coot {
namespace ideal_rna {

void delete_o2_prime_gemmi(gemmi::Residue *res) {
  if (!res) return;

  bool deleted = false;
  auto it = std::remove_if(
    res->atoms.begin(),
    res->atoms.end(),
    [&deleted](const gemmi::Atom &atom) {
      std::string atname = atom.name;
      if (atname == "O2*" || atname == "O2'") {
        deleted = true;
        return true;
      }
      return false;
    }
  );

  if (deleted) {
    res->atoms.erase(it, res->atoms.end());
  }
}

}
}