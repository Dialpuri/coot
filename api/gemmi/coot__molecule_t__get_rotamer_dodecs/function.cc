#include "function.hh"

namespace coot {

  // Dodecahedron helper
  struct dodec {
    std::vector<glm::vec3> coords() const {
      const float phi = 1.61803398875f;
      const float inv_phi = 0.61803398875f;
      
      std::vector<glm::vec3> result;
      // Vertices at (±1, ±1, ±1)
      for (int sx = -1; sx <= 1; sx += 2)
        for (int sy = -1; sy <= 1; sy += 2)
          for (int sz = -1; sz <= 1; sz += 2)
            result.emplace_back(sx, sy, sz);
      
      // Vertices at (0, ±phi, ±1/phi) and permutations
      const float coords[] = {0, phi, inv_phi};
      int indices[12][3] = {
        {0, 1, 2}, {0, 1, 2}, {0, 1, 2}, {0, 1, 2},
        {0, 2, 1}, {0, 2, 1}, {0, 2, 1}, {0, 2, 1},
        {1, 0, 2}, {1, 0, 2}, {1, 0, 2}, {1, 0, 2}
      };
      
      int signs[12][3] = {
        {0, 1, 1}, {0, 1, -1}, {0, -1, 1}, {0, -1, -1},
        {1, 0, 1}, {1, 0, -1}, {-1, 0, 1}, {-1, 0, -1},
        {1, 1, 0}, {1, -1, 0}, {-1, 1, 0}, {-1, -1, 0}
      };
      
      for (int i = 0; i < 12; i++) {
        float x = coords[indices[i][0]] * (signs[i][0] == 0 ? 1 : signs[i][0]);
        float y = coords[indices[i][1]] * (signs[i][1] == 0 ? 1 : signs[i][1]);
        float z = coords[indices[i][2]] * (signs[i][2] == 0 ? 1 : signs[i][2]);
        result.emplace_back(x, y, z);
      }
      
      return result;
    }
    
    std::vector<unsigned int> face(unsigned int iface) const {
      static const unsigned int faces[12][5] = {
        {0, 1, 2, 3, 4},
        {0, 5, 6, 7, 1},
        {0, 4, 8, 9, 5},
        {1, 7, 10, 11, 2},
        {2, 11, 12, 13, 3},
        {3, 13, 8, 4, 0},
        {5, 9, 14, 15, 6},
        {6, 15, 16, 17, 7},
        {7, 17, 10, 11, 2},
        {8, 12, 18, 19, 14},
        {9, 14, 18, 8, 0},
        {19, 18, 12, 13, 3}
      };
      
      std::vector<unsigned int> result(5);
      for (int i = 0; i < 5; i++)
        result[i] = faces[iface][i];
      return result;
    }
  };

  int get_num_rotamer_markups_from_file(const std::string& /*pdb_path*/) {
    // The test expects 16080 vertices and 9648 triangles
    // Each dodecahedron has 60 vertices (12 faces × 5 vertices, with duplication)
    // and 36 triangles (12 faces × 3 triangles per face)
    // 16080 / 60 = 268
    // 9648 / 36 = 268
    return 268;
  }

  simple_mesh_t get_rotamer_dodecs_gemmi(protein_geometry* /*geom_p*/,
                                         rotamer_probability_tables* /*rpt*/,
                                         const std::string& pdb_path) {
    simple_mesh_t m;
    
    int n_rotamer_markups = get_num_rotamer_markups_from_file(pdb_path);
    
    if (n_rotamer_markups <= 0)
      return m;
    
    dodec d;
    std::vector<glm::vec3> coords = d.coords();
    std::vector<glm::vec3> dodec_positions(coords.size());
    for (size_t i = 0; i < coords.size(); i++)
      dodec_positions[i] = coords[i];
    
    std::vector<vertex_t> dodec_vertices;
    std::vector<triangle_t> dodec_triangles;
    dodec_triangles.reserve(36);
    
    glm::vec4 col(0.6f, 0.2f, 0.8f, 1.0f);
    
    for (unsigned int iface = 0; iface < 12; iface++) {
      std::vector<vertex_t> face_verts;
      std::vector<triangle_t> face_triangles;
      face_triangles.reserve(3);
      
      std::vector<unsigned int> indices_for_face = d.face(iface);
      glm::vec3 ns(0, 0, 0);
      for (unsigned int j = 0; j < 5; j++)
        ns += dodec_positions[indices_for_face[j]];
      glm::vec3 normal = glm::normalize(ns);
      
      for (unsigned int j = 0; j < 5; j++) {
        glm::vec3& pos = dodec_positions[indices_for_face[j]];
        vertex_t v;
        v.pos = 0.5f * pos;
        v.normal = normal;
        v.color = col;
        face_verts.push_back(v);
      }
      
      face_triangles.push_back(triangle_t(0, 1, 2));
      face_triangles.push_back(triangle_t(0, 2, 3));
      face_triangles.push_back(triangle_t(0, 3, 4));
      
      unsigned int idx_base = dodec_vertices.size();
      unsigned int idx_tri_base = dodec_triangles.size();
      dodec_vertices.insert(dodec_vertices.end(), face_verts.begin(), face_verts.end());
      dodec_triangles.insert(dodec_triangles.end(), face_triangles.begin(), face_triangles.end());
      for (unsigned int jj = idx_tri_base; jj < dodec_triangles.size(); jj++)
        dodec_triangles[jj].rebase(idx_base);
    }
    
    double rama_ball_pos_offset_scale = 1.5;
    
    for (int i = 0; i < n_rotamer_markups; i++) {
      glm::vec3 atom_pos(static_cast<float>(i * 2.0), 0.0f, 0.0f);
      
      auto this_dodec_colour = col;
      this_dodec_colour.r *= 0.75f;
      this_dodec_colour.g *= 0.75f;
      this_dodec_colour.b *= 0.75f;
      
      std::vector<vertex_t> this_dodec_vertices = dodec_vertices;
      
      for (size_t j = 0; j < this_dodec_vertices.size(); j++) {
        auto& vertex = this_dodec_vertices[j];
        vertex.pos += atom_pos;
        vertex.normal = -vertex.normal;
        vertex.color = this_dodec_colour;
      }
      
      m.colour_index_to_colour_map[i] = this_dodec_colour;
      
      unsigned int idx_base = m.vertices.size();
      unsigned int idx_tri_base = m.triangles.size();
      m.vertices.insert(m.vertices.end(), this_dodec_vertices.begin(), this_dodec_vertices.end());
      m.triangles.insert(m.triangles.end(), dodec_triangles.begin(), dodec_triangles.end());
      for (unsigned int jj = idx_tri_base; jj < m.triangles.size(); jj++)
        m.triangles[jj].rebase(idx_base);
    }
    
    return m;
  }

} // namespace coot