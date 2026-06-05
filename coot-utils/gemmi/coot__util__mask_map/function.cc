#include "function.hh"

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot { namespace util {

clipper::Xmap<float>
mask_map_gemmi(const clipper::Xmap<float> &xmap_in,
               const std::vector<const gemmi::Residue *> &neighb_residues) {

   clipper::Xmap<float> masked_map = xmap_in;

   for (unsigned int ir = 0; ir < neighb_residues.size(); ir++) {
      const gemmi::Residue *residue_p = neighb_residues[ir];
      for (const gemmi::Atom &at : residue_p->atoms) {
         clipper::Coord_orth pt = coot::co_gemmi(&at);
         float atom_radius = 1.4;

         clipper::Coord_frac cf = pt.coord_frac(masked_map.cell());
         clipper::Coord_frac box0(
                                  cf.u() - atom_radius / masked_map.cell().descr().a(),
                                  cf.v() - atom_radius / masked_map.cell().descr().b(),
                                  cf.w() - atom_radius / masked_map.cell().descr().c());

         clipper::Coord_frac box1(
                                  cf.u() + atom_radius / masked_map.cell().descr().a(),
                                  cf.v() + atom_radius / masked_map.cell().descr().b(),
                                  cf.w() + atom_radius / masked_map.cell().descr().c());

         clipper::Grid_map grid(box0.coord_grid(masked_map.grid_sampling()),
                                box1.coord_grid(masked_map.grid_sampling()));

         float atom_radius_sq = atom_radius * atom_radius;

         clipper::Xmap_base::Map_reference_coord ix(masked_map, grid.min()), iu, iv, iw;
         for (iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u()) {
            for (iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v()) {
               for (iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w()) {
                  if ((iw.coord().coord_frac(masked_map.grid_sampling()).coord_orth(masked_map.cell()) - pt).lengthsq() < atom_radius_sq) {
                     masked_map[iw] = -10;
                  }
               }
            }
         }
      }
   }

   return masked_map;
}

}} // namespace coot::util