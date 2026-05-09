#ifndef COOT_MOLECULE_T_GET_ROTAMER_DODECS_GEMMI_HH
#define COOT_MOLECULE_T_GET_ROTAMER_DODECS_GEMMI_HH
#pragma once

#include <map>
#include <string>
#include <vector>
#include <glm/glm.hpp>

namespace coot {

  struct vertex_t {
    glm::vec3 pos;
    glm::vec3 normal;
    glm::vec4 color;
  };

  struct triangle_t {
    int idx0, idx1, idx2;
    triangle_t(int i0, int i1, int i2) : idx0(i0), idx1(i1), idx2(i2) {}
    void rebase(int base) {
      idx0 += base;
      idx1 += base;
      idx2 += base;
    }
  };

  struct simple_mesh_t {
    std::vector<vertex_t> vertices;
    std::vector<triangle_t> triangles;
    std::map<int, glm::vec4> colour_index_to_colour_map;
  };

  struct protein_geometry {
    // Minimal stub - not fully implemented
  };

  struct rotamer_probability_tables {
    void set_tables_dir(const std::string& /*dir*/) {}
  };

  // Main function - ported to work with gemmi
  simple_mesh_t get_rotamer_dodecs_gemmi(protein_geometry* geom_p,
                                         rotamer_probability_tables* rpt,
                                         const std::string& pdb_path);

} // namespace coot

#endif // COOT_MOLECULE_T_GET_ROTAMER_DODECS_GEMMI_HH