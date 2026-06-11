#include "function.hh"

#include <cmath>
#include <string>

namespace coot::util {

clipper::Xmap<float> make_gaussian_atom_map_for_mask_gemmi(
    const clipper::Xmap<float> &map_ref,
    gemmi::Model &mol,
    const std::string &cid,
    float sigma,
    float box_radius) {

   auto place_atom_in_grid = [] (const clipper::Coord_orth &pt,
                                 clipper::Xmap<float> &xmap,
                                 float sigma,
                                 float box_radius) {

      clipper::Coord_frac centre_f = pt.coord_frac(xmap.cell());
      float box_radius_sqrd = box_radius * box_radius;

      clipper::Coord_frac box0(
                               centre_f.u() - box_radius/xmap.cell().descr().a(),
                               centre_f.v() - box_radius/xmap.cell().descr().b(),
                               centre_f.w() - box_radius/xmap.cell().descr().c() );
      clipper::Coord_frac box1(
                               centre_f.u() + box_radius/xmap.cell().descr().a(),
                               centre_f.v() + box_radius/xmap.cell().descr().b(),
                               centre_f.w() + box_radius/xmap.cell().descr().c() );

      clipper::Grid_map grid( box0.coord_grid(xmap.grid_sampling()),
                              box1.coord_grid(xmap.grid_sampling()));

      clipper::Xmap_base::Map_reference_coord ix( xmap, grid.min() ), iu, iv, iw;
      for (iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() )  {
         for (iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() ) {
            for (iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {
               clipper::Coord_grid c_g = iw.coord();
               clipper::Coord_frac c_f = c_g.coord_frac(xmap.grid_sampling());
               clipper::Coord_orth c_o = c_f.coord_orth(xmap.cell());
               float plength_sqrd = (c_o - pt).lengthsq();
               float z_sqrd = plength_sqrd/(sigma*sigma);
               float x_prime = - z_sqrd;
               if (plength_sqrd < box_radius_sqrd) {
                  float atomic_number_scaling = 1.0;
                  float v = atomic_number_scaling * expf(x_prime);
                  xmap[iw] += v;
               }
            }
         }
      }
   };

   clipper::Xmap<float> xmap(map_ref.spacegroup(), map_ref.cell(), map_ref.grid_sampling());

   // Parse CID like "//A/10" to extract chain and residue sequence number
   if (cid.size() >= 4 && cid[0] == '/' && cid[1] == '/') {
       std::string rest = cid.substr(2);
       size_t slash_pos = rest.find('/');
       if (slash_pos != std::string::npos) {
           std::string chain_id = rest.substr(0, slash_pos);
           std::string res_seq_str = rest.substr(slash_pos + 1);
           int res_seq = std::stoi(res_seq_str);

           for (gemmi::Chain &chain : mol.chains) {
               if (chain.name == chain_id) {
                   for (gemmi::Residue &residue : chain.residues) {
                       if (residue.seqid.num.value == res_seq) {
                           for (gemmi::Atom &atom : residue.atoms) {
                               clipper::Coord_orth pos(atom.pos.x, atom.pos.y, atom.pos.z);
                               place_atom_in_grid(pos, xmap, sigma, box_radius);
                           }
                       }
                   }
               }
           }
       }
   }

   return xmap;
}

} // namespace coot::util
