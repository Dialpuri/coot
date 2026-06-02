#pragma once

#include <gemmi/model.hpp>
#include <gemmi/math.hpp>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <vector>
#include <string>
#include "api/make-instanced-graphical-bonds.hh"
#include "api/instancing.hh"
#include "api/bond-colour.hh"
#include "coot-utils/gl-matrix.h"
#include "coot-utils/oct.hh"
#include "coot-utils/cylinder.hh"
#include "coot-utils/ortep.hh"
#include "coot/analysis/chi-squared.hh"
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
    bool render_as_aniso;
};

// gemmi-compatible container
class graphical_bonds_container_gemmi {
public:
    int n_consolidated_atom_centres;
    graphical_bonds_points_list<graphical_bonds_atom_info_t_gemmi> *consolidated_atom_centres;
};

void make_instanced_graphical_bonds_spherical_atoms_gemmi(
    coot::instanced_mesh_t &m,
    const graphical_bonds_container_gemmi &gbc,
    coot::api_bond_colour_t bonds_box_type,
    float base_atom_radius,
    float base_bond_radius,
    bool render_atoms_as_aniso,
    float aniso_probability,
    bool render_aniso_atoms_as_ortep,
    bool render_aniso_atoms_as_empty,
    unsigned int num_subdivisions,
    const std::vector<glm::vec4> &colour_table);