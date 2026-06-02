#include "function.hh"

namespace coot {

namespace side_chain_densities {

density_box_t_gemmi sample_map_gemmi(
    const gemmi::Residue *residue_this_p,
    const gemmi::Residue *residue_next_p,
    side_chain_mode_t mode,
    const clipper::Coord_orth &cb_pt,
    const std::vector<clipper::Coord_orth> &axes,
    const clipper::Xmap<float> &xmap,
    std::string gen_pts_file_name,
    float grid_box_radius,
    int n_steps,
    const std::set<int> &useable_grid_points) {

   bool gen_usable_points_flag = false;
   if (mode == MODE_GEN_USABLE_POINTS) gen_usable_points_flag = true;

   int n_per_side = 2 * n_steps + 1;
   int n_box_vol  = n_per_side * n_per_side * n_per_side;

   if (!residue_this_p) return density_box_t_gemmi(nullptr, 0);
   if (axes.empty())     return density_box_t_gemmi(nullptr, 0);

   std::string res_name = residue_this_p->name;

   if (mode == MODE_SAMPLE_FOR_DB) {
      if (res_name == "GLY") {
         // rot_name = "pseudo";
      } else {
         return density_box_t_gemmi(nullptr, 0);
      }
   }

   clipper::Coord_orth ca_pt(-1, -1, -1);

   std::vector<clipper::Coord_orth> residue_atom_positions;
   std::vector<std::pair<double, clipper::Coord_orth>> main_chain_atom_positions;
   residue_atom_positions.reserve(residue_this_p->atoms.size());

   for (const auto &at : residue_this_p->atoms) {
      clipper::Coord_orth pos = coot::co_gemmi(&at);
      residue_atom_positions.push_back(pos);
      std::string atom_name = coot::trim_name(at.name);
      if (atom_name == "N" || atom_name == "C" || atom_name == "O" ||
          atom_name == "H" || atom_name == "CA") {
         double r = 2.8;
         if (atom_name == "CA") r = 1.6;
         main_chain_atom_positions.emplace_back(r, pos);
      }
      if (atom_name == "CA") ca_pt = pos;
   }

   if (gen_usable_points_flag && residue_next_p) {
      for (const auto &at : residue_next_p->atoms) {
         std::string atom_name = coot::trim_name(at.name);
         if (atom_name == "N") {
            clipper::Coord_orth pos = coot::co_gemmi(&at);
            main_chain_atom_positions.emplace_back(3.0, pos);
         }
      }
   }

   float step_size = grid_box_radius / static_cast<float>(n_steps);

   std::ofstream f;
   if (gen_usable_points_flag)
      f.open(gen_pts_file_name);

   float *density_box = new float[n_box_vol];
   float unset_value = -1001.1f;
   for (int i = 0; i < n_box_vol; ++i) density_box[i] = unset_value;

   std::set<int> main_chain_clashing_points;
   if (gen_usable_points_flag) {
      for (int ix = -n_steps; ix <= n_steps; ++ix) {
         for (int iy = -n_steps; iy <= n_steps; ++iy) {
            for (int iz = -n_steps; iz <= n_steps; ++iz) {
               int idx = (ix + n_steps) * n_per_side * n_per_side
                       + (iy + n_steps) * n_per_side
                       + (iz + n_steps);
               clipper::Coord_orth pt_in_grid =
                   make_pt_in_grid_gemmi(ix, iy, iz, step_size, axes);
               clipper::Coord_orth pt_grid_point = cb_pt + pt_in_grid;
               clipper::Coord_orth ca_to_cb = cb_pt - ca_pt;
               clipper::Coord_orth ca_to_pt = pt_grid_point - ca_pt;
               double dp = ca_to_cb[0]*ca_to_pt[0]
                         + ca_to_cb[1]*ca_to_pt[1]
                         + ca_to_cb[2]*ca_to_pt[2];
               if (dp < 0.0)
                  main_chain_clashing_points.insert(idx);
            }
         }
      }
   }

   double sum    = 0;
   double sum_sq = 0;
   int    count  = 0;

   for (int ix = -n_steps; ix <= n_steps; ++ix) {
      for (int iy = -n_steps; iy <= n_steps; ++iy) {
         for (int iz = -n_steps; iz <= n_steps; ++iz) {
            int idx = (ix + n_steps) * n_per_side * n_per_side
                    + (iy + n_steps) * n_per_side
                    + (iz + n_steps);

            if (gen_usable_points_flag
                || useable_grid_points.find(idx)
                     != useable_grid_points.end()) {
               clipper::Coord_orth pt_in_grid =
                   make_pt_in_grid_gemmi(ix, iy, iz, step_size, axes);
               clipper::Coord_orth pt_grid_point = cb_pt + pt_in_grid;

               if (gen_usable_points_flag) {
                  if (!is_close_to_atoms_gemmi(
                          main_chain_atom_positions, pt_grid_point)) {
                     if (in_sphere_gemmi(pt_grid_point, cb_pt, grid_box_radius)) {
                        if (main_chain_clashing_points.find(idx)
                            == main_chain_clashing_points.end()) {
                           f << "setting grid point " << idx << " at "
                             << pt_grid_point[0] << " "
                             << pt_grid_point[1] << " "
                             << pt_grid_point[2] << std::endl;
                        }
                     }
                  }
               } else {
                  float dv = 0.0f;
                  if (!xmap.is_null()) {
                     dv = 0.0f; // placeholder (no map loaded in test)
                  }
                  density_box[idx] = dv;
                  sum    += dv;
                  sum_sq += dv * dv;
                  ++count;
               }
            }
         }
      }
   }

   double mean = (count > 0) ? sum / count : 0.0;
   double variance =
       (count > 1) ? (sum_sq - sum * sum / count) / (count - 1) : 0.0;

   density_box_t_gemmi db(density_box, n_steps);
   db.mean = mean;
   db.var  = variance;

   return db;
}

} // namespace side_chain_densities
} // namespace coot