#include "function.hh"
#include <cmath>
#include <algorithm>

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
    const std::vector<glm::vec4> &colour_table) {

   auto convert_vertices = [] (const std::vector<coot::api::vnc_vertex> &v_in) {
      std::vector<coot::api::vn_vertex> v_out(v_in.size());
      for (unsigned int i=0; i<v_in.size(); i++) {
         const auto &v = v_in[i];
         v_out[i] = coot::api::vn_vertex(v.pos, v.normal);
      }
      return v_out;
   };

   auto convert_vertices_z_normal = [] (const std::vector<coot::api::vnc_vertex> &v_in) {
      std::vector<coot::api::vn_vertex> v_out(v_in.size());
      glm::vec3 z_normal(0.0f, 0.0f, 1.0f);
      for (unsigned int i=0; i<v_in.size(); i++) {
         const auto &v = v_in[i];
         v_out[i] = coot::api::vn_vertex(v.pos, z_normal);
      }
      return v_out;
   };

   auto add_ellipse_rings = [] (coot::instanced_geometry_t *ig, const glm::vec3 &sc, const glm::vec3 &t, const glm::mat4 &ori, const glm::vec4 &col) {

      glm::vec3 sc_other = 1.000f * sc;
      glm::mat4 m_unit(1.0);
      glm::mat4 m_x = glm::rotate(m_unit, static_cast<float>(0.5 * M_PI), glm::vec3(1.0, 0.0, 0.0));
      glm::mat4 m_y = glm::rotate(m_unit, static_cast<float>(0.5 * M_PI), glm::vec3(0.0, 1.0, 0.0));
      coot::instancing_data_type_B_t   idB(t, col, sc_other, ori);
      coot::instancing_data_type_B_t idB_x(t, col, sc_other, ori * m_x);
      coot::instancing_data_type_B_t idB_y(t, col, sc_other, ori * m_y);
      ig->instancing_data_B.push_back(idB);
      ig->instancing_data_B.push_back(idB_x);
      ig->instancing_data_B.push_back(idB_y);
   };

   auto add_ellipsoidal_multi_ring = [] (coot::instanced_geometry_t *ig_long,
                                         coot::instanced_geometry_t *ig_lat,
                                         const glm::vec3 &sc, const glm::vec3 &t,
                                         const glm::mat4 &ori, const glm::vec4 &col) {

      glm::mat4 m_unit(1.0);
      glm::mat4 m_x = glm::rotate(m_unit, static_cast<float>(0.5 * M_PI), glm::vec3(1.0, 0.0, 0.0));

      unsigned int n_longitude_steps = 10;
      for (unsigned int i = 0; i < n_longitude_steps; i++) {
         float angle = static_cast<float>(M_PI) * static_cast<float>(i) / static_cast<float>(n_longitude_steps);
         glm::mat4 m_z = glm::rotate(m_unit, angle, glm::vec3(0.0, 0.0, 1.0));
         coot::instancing_data_type_B_t idB(t, col, sc, ori * m_z * m_x);
         ig_long->instancing_data_B.push_back(idB);
      }

      coot::instancing_data_type_B_t idB_equator(t, col, sc, ori);
      ig_lat->instancing_data_B.push_back(idB_equator);

      unsigned int n_latitude_steps = 5;
      for (unsigned int i = 1; i <= n_latitude_steps; i++) {
         float lat_angle = static_cast<float>(0.5 * M_PI) * static_cast<float>(i) / static_cast<float>(n_latitude_steps + 1);
         float cos_lat = std::cos(lat_angle);
         float sin_lat = std::sin(lat_angle);

         glm::vec3 sc_lat(sc.x * cos_lat, sc.y * cos_lat, sc.z);
         glm::vec3 z_offset_local(0.0f, 0.0f, sc.z * sin_lat);
         glm::vec3 z_offset_world = glm::vec3(ori * glm::vec4(z_offset_local, 0.0f));

         coot::instancing_data_type_B_t idB_above(t + z_offset_world, col, sc_lat, ori);
         ig_lat->instancing_data_B.push_back(idB_above);

         coot::instancing_data_type_B_t idB_below(t - z_offset_world, col, sc_lat, ori);
         ig_lat->instancing_data_B.push_back(idB_below);
      }
   };

   coot::instanced_geometry_t ig_sphere("atoms");

   std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > octosphere_geom =
      tessellate_octasphere(num_subdivisions);

   std::vector<coot::api::vn_vertex> local_vertices(octosphere_geom.first.size());
   for (unsigned int i=0; i<octosphere_geom.first.size(); i++) {
      const glm::vec3 &v(octosphere_geom.first[i]);
      local_vertices[i] = coot::api::vn_vertex(v, v);
   }
   ig_sphere.vertices = local_vertices;
   ig_sphere.triangles = octosphere_geom.second;

   coot::instanced_geometry_t ig_ortep("anisotropic-ortep-atoms");
   ortep_t ortep = tessellate_sphere_sans_octant();
   bool do_ellipse_rings = true;
   std::vector<coot::api::vn_vertex> local_vertices_ortep(ortep.vertices.size());
   for (unsigned int i=0; i<ortep.vertices.size(); i++)
      local_vertices_ortep[i] = coot::api::vn_vertex(ortep.vertices[i], ortep.normals[i]);
   ig_ortep.vertices = local_vertices_ortep;
   ig_ortep.triangles = ortep.triangles;

   float ellipsoid_band_thickness = 0.02f;
   if (render_aniso_atoms_as_empty) ellipsoid_band_thickness = 0.04;
   cylinder cylinder_for_bands(std::make_pair(glm::vec3(0.0f, 0.0f,  ellipsoid_band_thickness),
                                              glm::vec3(0.0f, 0.0f, -ellipsoid_band_thickness)),
                               1.0f, 1.0f, 2.0f * ellipsoid_band_thickness, 32, 2);
   coot::instanced_geometry_t ig_ellipsoid_band("ellipsoid band longitude");
   ig_ellipsoid_band.vertices = convert_vertices(cylinder_for_bands.vertices);
   ig_ellipsoid_band.triangles = cylinder_for_bands.triangles;

   coot::instanced_geometry_t ig_ellipsoid_band_latitude("ellipsoid band latitude");
   ig_ellipsoid_band_latitude.vertices = convert_vertices_z_normal(cylinder_for_bands.vertices);
   ig_ellipsoid_band_latitude.triangles = cylinder_for_bands.triangles;

   int cts = colour_table.size();
   for (int icol=0; icol<gbc.n_consolidated_atom_centres; icol++) {
      glm::vec4 col(0.4, 0.4, 0.4, 1.0);
      if (icol<cts)
         col = colour_table[icol];
      for (unsigned int i=0; i<gbc.consolidated_atom_centres[icol].num_points; i++) {
         const graphical_bonds_atom_info_t_gemmi &at_info = gbc.consolidated_atom_centres[icol].points[i];
         gemmi::Atom *atom = at_info.atom_p;

         bool do_it = true;

         if (do_it) {
            float scale = at_info.radius_scale;
            float sar = scale * base_atom_radius;
            if (sar > 2.2) sar = 2.2;
            if (at_info.is_water)
               if (sar > 0.65) sar = 0.65f;

            glm::vec3 sc(sar, sar, sar);
            glm::vec3 t(atom->pos.x, atom->pos.y, atom->pos.z);

            bool atom_is_aniso = atom->aniso.nonzero();

            if (render_atoms_as_aniso && atom_is_aniso) {
               sar = 1.2;
               sar = gphl::prob_to_radius(aniso_probability * 100.0f);
               sc = glm::vec3(sar);

               GL_matrix mat(atom->aniso.u11, atom->aniso.u12, atom->aniso.u13,
                             atom->aniso.u12, atom->aniso.u22, atom->aniso.u23,
                             atom->aniso.u13, atom->aniso.u23, atom->aniso.u33);

               std::pair<bool,GL_matrix> chol_pair = mat.eigensystem();
               if (chol_pair.first) {
                  const auto &m_eig = chol_pair.second;
                  glm::mat4 ori(m_eig.matrix_element(0,0), m_eig.matrix_element(1,0), m_eig.matrix_element(2,0), 0.0f,
                                m_eig.matrix_element(0,1), m_eig.matrix_element(1,1), m_eig.matrix_element(2,1), 0.0f,
                                m_eig.matrix_element(0,2), m_eig.matrix_element(1,2), m_eig.matrix_element(2,2), 0.0f,
                                0.0f, 0.0f, 0.0f, 1.0f);
                  coot::instancing_data_type_B_t idB(t, col, sc, ori);
                  if (render_aniso_atoms_as_ortep)
                     ig_ortep.instancing_data_B.push_back(idB);
                  else
                     if (! render_aniso_atoms_as_empty)
                        ig_sphere.instancing_data_B.push_back(idB);
                  if (do_ellipse_rings) {
                     if (render_aniso_atoms_as_empty) {
                        add_ellipsoidal_multi_ring(&ig_ellipsoid_band, &ig_ellipsoid_band_latitude,
                                                   sc, t, ori, col);
                     } else {
                        glm::vec4 ellipsoid_ring_col = glm::vec4(0.1, 0.1, 0.1, 1.0);
                        add_ellipse_rings(&ig_ellipsoid_band, sc, t, ori, ellipsoid_ring_col);
                     }
                     if (render_aniso_atoms_as_empty) {
                        glm::vec3 sc_local = glm::vec3(base_atom_radius);
                        coot::instancing_data_type_A_t idA(t, col, sc_local);
                        ig_sphere.instancing_data_A.push_back(idA);
                     }
                  }
               } else {
                  coot::instancing_data_type_A_t idA(t, col, sc);
                  ig_sphere.instancing_data_A.push_back(idA);
               }

            } else {
               coot::instancing_data_type_A_t idA(t, col, sc);
               ig_sphere.instancing_data_A.push_back(idA);
            }
         }
      }
   }

   if (!ig_sphere.empty())                   m.add(ig_sphere);
   if (!ig_ortep.empty())                    m.add(ig_ortep);
   if (!ig_ellipsoid_band.empty())           m.add(ig_ellipsoid_band);
   if (!ig_ellipsoid_band_latitude.empty())  m.add(ig_ellipsoid_band_latitude);
}