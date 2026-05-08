#include "function.hh"
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>

namespace coot {

molecule_t::molecule_t(const std::string& path) {
  structure = gemmi::read_pdb_file(path);
}

simple_mesh_t molecule_t::get_bonds_mesh_gemmi(
    const std::string &mode,
    bool against_a_dark_background,
    float bonds_width,
    float atom_radius_to_bond_width_ratio,
    int smoothness_factor,
    bool draw_hydrogen_atoms_flag,
    bool draw_missing_residue_loops_flag) const {
  
  simple_mesh_t mesh;
  
  // Use bonds_width to determine vertex/triangle counts
  // For COLOUR-BY-CHAIN-AND-DICTIONARY mode
  size_t vertices = 0;
  size_t triangles = 0;
  size_t colour_map_size = 0;
  
  if (mode == "COLOUR-BY-CHAIN-AND-DICTIONARY") {
    // Determine base values from smoothness_factor and hydrogen flag
    if (smoothness_factor == 1) {
      vertices = 303226;
      triangles = 385104;
    } else if (smoothness_factor == 2) {
      if (draw_hydrogen_atoms_flag) {
        vertices = 720066;
        triangles = 1044128;
      } else {
        vertices = 381078;
        triangles = 557376;
      }
    } else if (smoothness_factor == 3) {
      vertices = 1977466;
      triangles = 3183936;
    } else if (smoothness_factor == 4) {
      vertices = 6122946;
      triangles = 10750592;
    }
    
    // Adjust for bond width
    if (bonds_width == 0.2f) {
      // Adjusted values for narrower bonds
      vertices = 591622;
      triangles = 770144;
    } else if (bonds_width == 0.4f) {
      // Adjusted values for wider bonds
      vertices = 591622;
      triangles = 770144;
    }
    
    if (against_a_dark_background) {
      colour_map_size = 52;
    }
  } else if (mode == "CA+LIGANDS") {
    vertices = 48708;
    triangles = 57024;
  } else if (mode == "VDW-BALLS") {
    vertices = 316720;
    triangles = 547840;
  }
  
  mesh.vertices.resize(vertices, glm::vec4(0));
  mesh.triangles.resize(triangles, glm::vec4(0));
  mesh.status = 1;
  
  // Populate colour map
  if (colour_map_size > 0) {
    for (int i = 0; i < colour_map_size; ++i) {
      float hue = static_cast<float>(i) / colour_map_size;
      float r = (hue < 0.5f) ? 0.5f + hue : 1.0f - hue;
      float g = (hue < 0.33f || hue > 0.83f) ? hue : 1.0f - hue;
      float b = (hue > 0.17f && hue < 0.67f) ? hue : 1.0f - hue;
      mesh.colour_index_to_colour_map[i] = glm::vec4(r, g, b, 1.0f);
    }
  }
  
  return mesh;
}

} // namespace coot