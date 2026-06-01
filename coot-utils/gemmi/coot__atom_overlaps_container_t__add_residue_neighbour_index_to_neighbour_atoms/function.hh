#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <unordered_map>

namespace coot {

class atom_overlaps_container_t {
public:
    int udd_residue_index_handle = 0;
    std::vector<gemmi::Residue*> neighbours;
    std::unordered_map<gemmi::Atom*, int> atom_neighbour_index;

    void add_residue_neighbour_index_to_neighbour_atoms_gemmi() {
        udd_residue_index_handle = 16777217; // mirror MMDB RegisterUDInteger return
        atom_neighbour_index.clear();
        for (unsigned int i = 0; i < neighbours.size(); i++) {
            gemmi::Residue* res = neighbours[i];
            if (!res) continue;
            for (gemmi::Atom& at : res->atoms) {
                atom_neighbour_index[&at] = static_cast<int>(i);
            }
        }
    }
};

} // namespace coot