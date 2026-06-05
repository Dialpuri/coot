#include "function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"
#include <clipper/contrib/skeleton.h>
#include <cmath>
#include <queue>
#include <vector>

clipper::Xmap<float>
coot::util::make_map_mask_gemmi(const clipper::Spacegroup &space_group,
                                const clipper::Cell &cell,
                                const clipper::Grid_sampling &grid_sampling,
                                std::vector<const gemmi::Atom*> atoms,
                                float radius,
                                float smooth) {

   float outer_space_f = -1.1f;
   clipper::Xmap<float> xmap(space_group, cell, grid_sampling);
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = xmap.first(); !ix.last(); ix.next())
      xmap[ix] = outer_space_f;

   for (size_t iat = 0; iat < atoms.size(); iat++) {
      const gemmi::Atom* at = atoms[iat];

      float atom_radius_sq = radius * radius;
      clipper::Coord_orth pt = coot::co_gemmi(at);
      clipper::Coord_frac cf = pt.coord_frac(cell);
      clipper::Coord_frac box0(
                               cf.u() - radius/cell.descr().a(),
                               cf.v() - radius/cell.descr().b(),
                               cf.w() - radius/cell.descr().c());
      clipper::Coord_frac box1(
                               cf.u() + radius/cell.descr().a(),
                               cf.v() + radius/cell.descr().b(),
                               cf.w() + radius/cell.descr().c());
      clipper::Grid_map grid(box0.coord_grid(grid_sampling),
                             box1.coord_grid(grid_sampling));

      if (smooth > 0.0) {
         clipper::Xmap_base::Map_reference_coord ix(xmap, grid.min() ), iu, iv, iw;
         for (iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() ) {
            for ( iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() ) {
               for ( iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {
                  clipper::Coord_orth gp_orth = iw.coord().coord_frac(grid_sampling).coord_orth(cell);
                  float dd = (gp_orth - pt).lengthsq();
                  const float &current_v = xmap[iw];
                  if (current_v < 1.0) {
                     if (dd < atom_radius_sq) {
                        float d = sqrtf(dd);
                        float v_delta = 1.0 - d/radius;
                        xmap[iw] += v_delta;
                        if (xmap[iw] > 1.0) xmap[iw] = 1.0;
                     }
                  }
               }
            }
         }

      } else {
         // no smoothing
         clipper::Xmap_base::Map_reference_coord ix(xmap, grid.min() ), iu, iv, iw;
         for (iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() ) {
            for ( iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() ) {
               for ( iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {
                  float dd = (iw.coord().coord_frac(grid_sampling).coord_orth(cell) - pt).lengthsq();
                  if (dd < atom_radius_sq) {
                     xmap[iw] = 1.0;
                  }
               }
            }
         }
      }
   }

   // remove voids
   bool continue_removing_voids = true;
   int count = 0;
   float cut_off = 1.0f;

   clipper::Xmap<short int> considered_map(space_group, cell, grid_sampling);
   clipper::Xmap<short int> outer_surface(space_group, cell, grid_sampling);
   for (ix = considered_map.first(); !ix.last(); ix.next()) considered_map[ix] = 0;
   for (ix = outer_surface.first();  !ix.last(); ix.next())  outer_surface[ix] = 0;

   clipper::Skeleton_basic::Neighbours neighb(xmap);

   while (continue_removing_voids) {

      continue_removing_voids = false;

      std::queue<clipper::Coord_grid> q;
      for (ix = xmap.first(); !ix.last(); ix.next()) {
         const float &xmix = xmap[ix];
         if (xmix < cut_off) {
            if (xmix > 0.0f) {
               q.push(ix.coord());
            }
         }

         if (! q.empty()) {
            if (considered_map.get_data(ix.coord()) == 1) {
            } else {
               considered_map[ix] = 1;
               std::vector<clipper::Coord_grid> cluster_grid_points;
               while (! q.empty()) {
                  clipper::Coord_grid c_g_start = q.front();
                  q.pop();
                  cluster_grid_points.push_back(c_g_start);
                  for (int i=0; i<neighb.size(); i++) {
                     clipper::Coord_grid c_g = c_g_start + neighb[i];
                     if (considered_map.get_data(c_g) == 0) {
                        float f = xmap.get_data(c_g);
                        if (f < cut_off) {
                           considered_map.set_data(c_g, 1);
                           cluster_grid_points.push_back(c_g);
                           if (f > 0.0f)
                              q.push(c_g);
                        }
                     }
                  }
               }

               if (cluster_grid_points.size() < 100) {
                  for (unsigned int ii=0; ii<cluster_grid_points.size(); ii++) {
                     const auto &gp = cluster_grid_points[ii];
                     xmap.set_data(gp, 1.0f);
                  }
                  continue_removing_voids = true;
               } else {
                  for (unsigned int ii=0; ii<cluster_grid_points.size(); ii++) {
                     const auto &gp = cluster_grid_points[ii];
                     outer_surface.set_data(gp, 1);
                  }
               }
            }
         }
      }
      count += 1;
   }

   // what's left is voids with values of actually zero
   continue_removing_voids = true;
   while (continue_removing_voids) {

      continue_removing_voids = false;

      for (ix = xmap.first(); !ix.last(); ix.next()) {
         const float &xmix = xmap[ix];
         if (xmix == 0.0f) {
            clipper::Coord_grid c_g_start = ix.coord();
            if (considered_map.get_data(c_g_start) == 1) {
            } else {
               bool has_non_zero_neighb = false;
               for (int i=0; i<neighb.size(); i++) {
                  clipper::Coord_grid c_g = c_g_start + neighb[i];
                  if (xmap.get_data(c_g) > 0.0) {
                     if (outer_surface.get_data(c_g) != 1) {
                        has_non_zero_neighb = true;
                        break;
                     }
                     if (has_non_zero_neighb) {
                        xmap.set_data(c_g, 1.0);
                        continue_removing_voids = true;
                     }
                  }
               }
            }
         }
      }
   }

   return xmap;
}