#pragma once

#include <gemmi/model.hpp>
#include <map>
#include <string>
#include <vector>

namespace coot::reduce {

struct atom_with_attached_Hs {
    enum hydrogen_t { METHYL, HYDROXYL, AMIDE, SULFHYDRYL, NONE };

    gemmi::Atom *atom;
    hydrogen_t type;
    std::vector<gemmi::Atom *> hydrogens;

    atom_with_attached_Hs(gemmi::Atom *at, hydrogen_t type_,
                          const std::vector<gemmi::Atom *> &attached_hydrogens)
        : atom(at), type(type_), hydrogens(attached_hydrogens) {}
};

class atoms_with_spinnable_Hs {
public:
    std::map<std::string, std::vector<atom_with_attached_Hs>> typed_atoms;

    void add_gemmi(gemmi::Atom *at, atom_with_attached_Hs::hydrogen_t type,
                   const std::vector<gemmi::Atom *> &attached_hydrogen_atoms) {
        std::string alt_loc;
        if (at && at->altloc != ' ' && at->altloc != '\0') {
            alt_loc = std::string(1, at->altloc);
        }
        atom_with_attached_Hs awaH(at, type, attached_hydrogen_atoms);
        typed_atoms[alt_loc].push_back(awaH);
    }
};

} // namespace coot::reduce
