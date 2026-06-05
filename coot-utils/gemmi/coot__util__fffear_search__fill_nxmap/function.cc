#include "function.hh"

#include <cmath>
#include <iostream>
#include <string>

#include <clipper/clipper.h>
#include <gemmi/model.hpp>

namespace coot {
namespace util {
namespace fffear_search {

int fill_nxmap_gemmi(clipper::NXmap<float>& nxmap,
                     const std::vector<const gemmi::Atom*>& atoms,
                     const clipper::Coord_orth& mid_point) {

   int n_atoms = static_cast<int>(atoms.size());

   clipper::ftype radius_ = 2.5;

   clipper::NXmap<float>::Map_reference_index im;
   nxmap = 0.0;

   int n_points = 0;
   clipper::Coord_grid g0, g1;
   g0 = clipper::Coord_map(nxmap.operator_orth_grid().rot() *
                           clipper::Vec3<>(radius_, radius_, radius_)).coord_grid();
   clipper::Grid_range gd(-g0, g0);
   clipper::Grid_range box(clipper::Coord_grid(0, 0, 0),
                           clipper::Coord_grid(nxmap.grid()) - clipper::Coord_grid(1, 1, 1));
   clipper::NXmap<float>::Map_reference_coord i0, iu, iv, iw;
   for (int i = 0; i < n_atoms; i++) {
      const gemmi::Atom* atom = atoms[i];
      if (atom) {
         clipper::Coord_orth p(atom->pos.x, atom->pos.y, atom->pos.z);
         p -= mid_point;
         clipper::AtomShapeFn sf(p, atom->element.name(),
                                 static_cast<double>(atom->b_iso),
                                 static_cast<double>(atom->occ));
         g0 = nxmap.coord_map(p).coord_grid() + gd.min();
         g1 = nxmap.coord_map(p).coord_grid() + gd.max();
         i0 = clipper::NXmap<float>::Map_reference_coord(nxmap, g0);
         float r;
         for (iu = i0; iu.coord().u() <= g1.u(); iu.next_u())
            for (iv = iu; iv.coord().v() <= g1.v(); iv.next_v())
               for (iw = iv; iw.coord().w() <= g1.w(); iw.next_w())
                  if (box.in_grid(iw.coord())) {
                     r = sf.rho(iw.coord_orth());
                     if (clipper::Util::isnan(r)) {
                        // caused by B factor of 0.00
                     } else {
                        nxmap[iw] += r;
                     }
                     n_points++;
                  }
      }
   }
   std::cout << "INFO:: Number of non-zero points in atom search map: " << n_points
             << std::endl;

   // debugging check for nans:
   clipper::NXmap<float>::Map_reference_index inx;
   float d;
   int n_nan = 0;
   int n_nx_points = 0;
   for (inx = nxmap.first(); !inx.last(); inx.next()) {
      n_nx_points++;
      d = nxmap[inx];
      if (clipper::Util::isnan(d)) {
         n_nan++;
      }
   }
   if (n_nan > 0) {
      std::cout << "----:: " << n_nan << " of " << n_nx_points
                << " map points were nans" << std::endl;
      std::cout << "----:: " << n_points << " were set to density values" << std::endl;
   }
   return n_points;
}

} // namespace fffear_search
} // namespace util
} // namespace coot