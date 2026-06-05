#pragma once

#include <gemmi/model.hpp>
#include <set>
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>

void add_fixed_atoms_from_non_bonded_neighbours_gemmi(
    const std::vector<gemmi::Residue *> &non_bonded_neighbour_residues,
    std::set<int> &fixed_atom_indices,
    const std::map<gemmi::Atom *, int> &atom_index_map)
{
    for (std::size_t jj = 0; jj < non_bonded_neighbour_residues.size(); jj++) {
        gemmi::Residue *residue_p = non_bonded_neighbour_residues[jj];
        for (gemmi::Atom &at : residue_p->atoms) {
            int idx = -1;
            if (atom_index_map.find(&at) == atom_index_map.end()) {
                std::cerr << "ERROR:: in add_fixed_atoms_from_non_bonded_neighbours_gemmi() "
                          << " bad atom_index_map for atom " << jj << std::endl;
            } else {
                idx = atom_index_map.find(&at)->second;
                if (std::find(fixed_atom_indices.begin(),
                              fixed_atom_indices.end(), idx) == fixed_atom_indices.end()) {
                    fixed_atom_indices.insert(idx);
                }
            }
        }
    }
}