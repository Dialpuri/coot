#pragma once

#include <gemmi/model.hpp>
#include <coot/api/instancing.hh>
#include <coot/coot-utils/vertex.hh>
#include <coot/coot-utils/oct.hh>
#include <coot/geometry/protein-geometry.hh>
#include <glm/glm.hpp>
#include <map>
#include <vector>

// Gemmi-compatible version of graphical_bonds_atom_info_t
// Uses gemmi::Atom* instead of mmdb::Atom*
class graphical_bonds_atom_info_t_gemmi {
public:
   bool is_hydrogen_atom;
   bool is_water;
   float radius_scale;
   glm::vec3 position;
   const gemmi::Atom *atom_p;
   int model_number;
   int atom_index;

   graphical_bonds_atom_info_t_gemmi() :
      is_hydrogen_atom(false), is_water(false), radius_scale(1.0f),
      atom_p(nullptr), model_number(-1), atom_index(-1) {}
};

// Gemmi-compatible version of graphical_bonds_consolidated_atom_centre
class graphical_bonds_consolidated_atom_centre_gemmi {
public:
   std::vector<graphical_bonds_atom_info_t_gemmi> points;
   unsigned int num_points;
   glm::vec4 colour;
   float radius;

   graphical_bonds_consolidated_atom_centre_gemmi() :
      num_points(0), colour(1.0f, 1.0f, 1.0f, 1.0f), radius(1.5f) {}
};

// Gemmi-compatible version of graphical_bonds_container
class graphical_bonds_container_gemmi {
public:
   int n_consolidated_atom_centres;
   std::vector<graphical_bonds_consolidated_atom_centre_gemmi> consolidated_atom_centres;

   graphical_bonds_container_gemmi() : n_consolidated_atom_centres(0) {}
};

namespace coot {

void make_graphical_bonds_spherical_atoms_with_vdw_radii_instanced_gemmi(
    coot::instanced_mesh_t &m,
    const graphical_bonds_container_gemmi &gbc,
    unsigned int num_subdivisions,
    const std::vector<glm::vec4> &colour_table,
    const coot::protein_geometry &geom,
    int imol_no) {

  coot::instanced_geometry_t ig("vdW Balls");

  // ----------------------- setup the vertices and triangles ---------------------

  std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > octosphere_geom =
     tessellate_octasphere(num_subdivisions);

  std::vector<coot::api::vn_vertex> local_vertices(octosphere_geom.first.size());
  for (unsigned int i=0; i<octosphere_geom.first.size(); i++) {
     const glm::vec3 &v(octosphere_geom.first[i]);
     local_vertices[i] = coot::api::vn_vertex(v, v);
  }
  ig.vertices = local_vertices;
  ig.triangles = octosphere_geom.second;

  // ----------------------- setup the instances ----------------------

  std::map<std::string, float> ele_to_radius_map;
  glm::mat4 unit(1.0);
  for (int icol=0; icol<gbc.n_consolidated_atom_centres; icol++) {
     glm::vec4 col = colour_table[icol];
     for (unsigned int i=0; i<gbc.consolidated_atom_centres[icol].num_points; i++) {
        const graphical_bonds_atom_info_t_gemmi &at_info = gbc.consolidated_atom_centres[icol].points[i];
        auto* at = at_info.atom_p;
        std::string ele(at->element.name());
        std::map<std::string, float>::const_iterator it = ele_to_radius_map.find(ele);
        float atom_radius = 1.0;
        if (it != ele_to_radius_map.end()) {
           atom_radius = it->second;
        } else {
           std::string atom_name(at->name);
           atom_radius = geom.get_vdw_radius(atom_name, "", imol_no, false);
           ele_to_radius_map[ele] = atom_radius;
        }

        glm::vec3 t(at->pos.x, at->pos.y, at->pos.z);
        glm::vec3 sc(atom_radius, atom_radius, atom_radius);
        coot::instancing_data_type_A_t id(t, col, sc);
        ig.instancing_data_A.push_back(id);
     }
  }
  m.add(ig);
}

} // namespace coot