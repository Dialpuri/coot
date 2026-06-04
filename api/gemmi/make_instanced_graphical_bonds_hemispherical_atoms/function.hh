#pragma once

#include <gemmi/model.hpp>
#include <glm/glm.hpp>
#include <vector>
#include "api/instancing.hh"
#include "api/bond-colour.hh"
#include "coords/graphical-bonds-container.hh"
#include "coot/coords/Cartesian.hh"

// gemmi-compatible atom info struct (atom_p uses gemmi::Atom* instead of mmdb::Atom*)
struct graphical_bonds_atom_info_t_gemmi {
    gemmi::Atom *atom_p;
    coot::Cartesian position;
    int atom_index;
    int model_number;
    bool is_hydrogen_atom;
    bool is_water;
    float radius_scale;
};

// gemmi-compatible container
class graphical_bonds_container_gemmi {
public:
    int n_consolidated_atom_centres;
    graphical_bonds_points_list<graphical_bonds_atom_info_t_gemmi> *consolidated_atom_centres;
};

void
make_instanced_graphical_bonds_hemispherical_atoms_gemmi(
    coot::instanced_mesh_t &m,
    const graphical_bonds_container_gemmi &gbc,
    coot::api_bond_colour_t bonds_box_type,
    float atom_radius,
    float bond_radius,
    unsigned int num_subdivisions,
    const std::vector<glm::vec4> &colour_table) {

   return; // 20230224-PE every atom is spherical for the moment.
}