#pragma once

#include <string>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>

namespace coot {

class restraints_container_t {
public:
    /// Gemmi port of get_N_index: returns the index of the "N" atom
    /// within the residue's atom vector, or -2 if not found.
    static int get_N_index_gemmi(const gemmi::Residue &residue) {
        int idx = -2;
        for (int i = 0; i < static_cast<int>(residue.atoms.size()); i++) {
            if (residue.atoms[i].name == "N") {
                idx = i;
                break;
            }
        }
        return idx;
    }
};

} // namespace coot