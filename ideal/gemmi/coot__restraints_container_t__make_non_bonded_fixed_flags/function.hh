#pragma once

#include <vector>
#include <set>
#include <algorithm>
#include <gemmi/model.hpp>

namespace coot {

struct restraints_container_t_gemmi {
    std::vector<gemmi::Atom> atom;
    std::vector<gemmi::Residue*> residue_for_atom;
    std::set<int> fixed_atom_indices;
    std::vector<gemmi::Residue*> non_bonded_neighbour_residues;
    size_t n_atoms = 0;

    std::vector<bool> make_non_bonded_fixed_flags_gemmi(int index1, int index2) const;
};

inline std::vector<bool>
restraints_container_t_gemmi::make_non_bonded_fixed_flags_gemmi(int index1, int index2) const {

    std::vector<bool> r(2, false);
    bool set_0 = false;
    bool set_1 = false;

    if (fixed_atom_indices.find(index1) != fixed_atom_indices.end()) {
        r[0] = true;
        set_0 = true;
    }
    if (fixed_atom_indices.find(index2) != fixed_atom_indices.end()) {
        r[1] = true;
        set_1 = true;
    }

    if (set_0 && set_1) {
        return r;
    }

    if (!set_0) {
        gemmi::Residue* res = residue_for_atom[index1];
        if (std::find(non_bonded_neighbour_residues.begin(),
                      non_bonded_neighbour_residues.end(),
                      res) != non_bonded_neighbour_residues.end())
            r[0] = true;
    }
    if (!set_1) {
        gemmi::Residue* res = residue_for_atom[index2];
        if (std::find(non_bonded_neighbour_residues.begin(),
                      non_bonded_neighbour_residues.end(),
                      res) != non_bonded_neighbour_residues.end())
            r[1] = true;
    }
    return r;
}

} // namespace coot