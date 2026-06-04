#pragma once

#include <vector>
#include <gemmi/model.hpp>

namespace coot {

struct restraints_container_t {
    std::vector<std::vector<int>> filtered_non_bonded_atom_indices;
    int n_atoms = 0;

    void filter_non_bonded_by_distance_gemmi(
        const std::vector<const gemmi::Atom*>& atoms,
        const std::vector<std::vector<int>>& non_bonded_atom_indices,
        double dist);
};

inline void coot::restraints_container_t::filter_non_bonded_by_distance_gemmi(
    const std::vector<const gemmi::Atom*>& atoms,
    const std::vector<std::vector<int>>& non_bonded_atom_indices,
    double dist)
{
    filtered_non_bonded_atom_indices.resize(non_bonded_atom_indices.size());

    double dist_lim2 = dist * dist;

    for (unsigned int i = 0; i < non_bonded_atom_indices.size(); i++) {
        for (unsigned int j = 0; j < non_bonded_atom_indices[i].size(); j++) {
            const gemmi::Atom& atom_1 = *atoms[i];
            const gemmi::Atom& atom_2 = *atoms[non_bonded_atom_indices[i][j]];

            gemmi::Vec3 diff = atom_1.pos - atom_2.pos;
            double dist2 = diff.length_sq();

            if (dist2 < dist_lim2) {
                // In MMDB, atom_2->GetUDData(udd_atom_index_handle, i_at_ind)
                // retrieves the atom's array index that was stored during setup.
                // This equals non_bonded_atom_indices[i][j].
                int i_at_ind = non_bonded_atom_indices[i][j];
                filtered_non_bonded_atom_indices[i].push_back(i_at_ind);
            }
        }
    }
}

} // namespace coot