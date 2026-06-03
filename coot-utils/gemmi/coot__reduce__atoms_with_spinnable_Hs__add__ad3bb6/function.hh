#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <map>

namespace coot {
namespace reduce {

struct atom_with_attached_Hs {
  enum class hydrogen_t {
    NONE = 0,
    HYDROXYL,
    SULFHYDRYL,
    AMIDE,
    METHYL
  };

  gemmi::Atom* atom;
  hydrogen_t type;
  std::vector<gemmi::Atom*> hydrogens;

  atom_with_attached_Hs(gemmi::Atom* at, hydrogen_t t, std::vector<gemmi::Atom*> h)
    : atom(at), type(t), hydrogens(std::move(h)) {}
};

struct atoms_with_spinnable_Hs {
  std::map<std::string, std::vector<atom_with_attached_Hs>> typed_atoms;

  void add_gemmi(gemmi::Atom* at,
                 atom_with_attached_Hs::hydrogen_t type,
                 gemmi::Atom* attached_hydrogen_atom) {
    std::string alt_loc(1, at->altloc);
    std::vector<gemmi::Atom*> v;
    v.push_back(attached_hydrogen_atom);
    atom_with_attached_Hs awaH(at, type, v);
    typed_atoms[alt_loc].push_back(awaH);
  }

  void add_gemmi(gemmi::Atom* at,
                 atom_with_attached_Hs::hydrogen_t type,
                 std::vector<gemmi::Atom*> attached_hydrogens) {
    std::string alt_loc(1, at->altloc);
    atom_with_attached_Hs awaH(at, type, attached_hydrogens);
    typed_atoms[alt_loc].push_back(awaH);
  }
};

} // namespace reduce
} // namespace coot