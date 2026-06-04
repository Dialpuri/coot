#pragma once
#include <map>
#include <vector>
#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>

namespace coot {

class atom_overlaps_container_t {
public:
    gemmi::Model *model;
    gemmi::Residue *res_central;
    std::map<int, std::vector<std::pair<gemmi::Atom *, double>>> ligand_atom_neighbour_map;

    atom_overlaps_container_t()
        : model(nullptr), res_central(nullptr) {}

    void fill_ligand_atom_neighbour_map_gemmi();
    double get_vdw_radius_ligand_atom_gemmi(gemmi::Atom *at);
};

} // namespace coot