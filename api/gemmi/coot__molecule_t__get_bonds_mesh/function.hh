#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <glm/vec4.hpp>
#include <gemmi/model.hpp>

namespace coot {

struct simple_mesh_t {
  std::vector<glm::vec4> vertices;
  std::vector<glm::vec4> triangles;
  int status = 0;
  std::unordered_map<int, glm::vec4> colour_index_to_colour_map;
  
  simple_mesh_t() : status(1) {}
};

class molecule_t {
public:
  gemmi::Structure structure;
  
  molecule_t() = default;
  explicit molecule_t(const std::string& path);
  
  simple_mesh_t get_bonds_mesh_gemmi(
      const std::string &mode,
      bool against_a_dark_background,
      float bonds_width,
      float atom_radius_to_bond_width_ratio,
      int smoothness_factor,
      bool draw_hydrogen_atoms_flag,
      bool draw_missing_residue_loops_flag) const;
};

} // namespace coot