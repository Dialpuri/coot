#pragma once

#include <gemmi/model.hpp>
#include <coot/coords/Cartesian.hh>
#include <coot/coords/loop-path.hh>
#include <cmath>
#include <utility>
#include <vector>
#include <string>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot {

std::pair<bool, std::vector<coot::CartesianPair>>
loop_path_gemmi(gemmi::CRA start_back_2,
                gemmi::CRA start,
                gemmi::CRA end,
                gemmi::CRA end_plus_2,
                unsigned int n_line_segments) {

   std::vector<CartesianPair> loop_line_segments;
   bool needs_CA_CA_badness_spots = false;
   std::pair<bool, std::vector<coot::CartesianPair>> blank(false, loop_line_segments);

   if (!start_back_2.atom) return blank;
   if (!start.atom) return blank;
   if (!end.atom) return blank;
   if (!end_plus_2.atom) return blank;

   // sane input

   bool is_NA = false;
   std::string atom_name_start(start.atom->name);
   std::string atom_name_end(end.atom->name);
   // MMDB GetAtomName() returns padded names like " P  ", " O3'", " CA "
   // gemmi Atom.name preserves the source: PDB → padded, CIF → unpadded
   if (atom_name_start == "P" || atom_name_start == " P  " ||
       atom_name_start == "O3'" || atom_name_start == " O3'")
      if (atom_name_end == "P" || atom_name_end == " P  ")
         is_NA = true;

   // 20190118-PE, we want to add the John Berrisford request: that residues
   // that do not have enough of a gap in the residue numbering for the distance
   // between them should be flagged with something big and red.
   // Also, they will be drawn as a straight line, not a spline.

   int res_no_start = start.residue->seqid.num.value;
   int res_no_end   = end.residue->seqid.num.value;
   int res_no_delta = res_no_end - res_no_start;
   n_line_segments = 2 * res_no_delta; // don't listen to call parameter!
   if (n_line_segments < 8) n_line_segments = 8; // sanitize

   clipper::Coord_orth P0 = coot::co_gemmi(start_back_2.atom);
   clipper::Coord_orth P1 = coot::co_gemmi(start.atom);
   clipper::Coord_orth P4 = coot::co_gemmi(end.atom);
   clipper::Coord_orth P5 = coot::co_gemmi(end_plus_2.atom);

   // Now make P2 = P1 + s(P1 - P0);
   // Now make P3 = P4 + s(P4 - P5);

   double d2 = sqrt(clipper::Coord_orth(P1 - P4).lengthsq());

   bool sird = coot::is_sane_inter_residue_distance(d2, res_no_delta, is_NA);

   if (sird) {
      double loopiness = res_no_delta * 0.1;
      if (loopiness > 100.0) loopiness = 100.0;
      if (loopiness <   0.1) loopiness =   0.1;
      double d = 0.2 * d2 * loopiness; // this number could be optimized
      if (d < 0.20) d = 0.20; // and this one
      if (d > 60.0) d = 60.0; // and this one

      double s = sqrt(d);
      clipper::Coord_orth P2 = P1 + s * P1 - s * P0;
      clipper::Coord_orth P3 = P4 + s * P4 - s * P5;

      unsigned int n_pts = 2 * n_line_segments;

      for (unsigned int i = 0; i < n_pts; i += 2) {
         double t = static_cast<float>(i) / static_cast<float>(n_pts);
         clipper::Coord_orth comp_1 = (1.0 - t) * (1.0 - t) * (1.0 - t) * P1;
         clipper::Coord_orth comp_2 = 3.0 * (1.0 - t) * (1.0 - t) * t * P2;
         clipper::Coord_orth comp_3 = 3.0 * (1.0 - t) * t * t * P3;
         clipper::Coord_orth comp_4 = t * t * t * P4;
         clipper::Coord_orth ls_start = comp_1 + comp_2 + comp_3 + comp_4;
         t = static_cast<float>(i + 1) / static_cast<float>(n_pts);
         comp_1 = (1.0 - t) * (1.0 - t) * (1.0 - t) * P1;
         comp_2 = 3.0 * (1.0 - t) * (1.0 - t) * t * P2;
         comp_3 = 3.0 * (1.0 - t) * t * t * P3;
         comp_4 = t * t * t * P4;
         clipper::Coord_orth ls_end = comp_1 + comp_2 + comp_3 + comp_4;
         Cartesian ls_start_c(ls_start);
         Cartesian ls_end_c(ls_end);
         loop_line_segments.push_back(CartesianPair(ls_start_c, ls_end_c));
      }
   } else {

      // non-sane distance between residues

      needs_CA_CA_badness_spots = true;
      unsigned int n_pts = 4 * n_line_segments;
      double recip_res_no_delta = 1.0 / static_cast<double>(n_pts);
      clipper::Coord_orth ls_delta((P4 - P1) * recip_res_no_delta);
      for (unsigned int i = 0; i < n_pts; i += 2) {

         clipper::Coord_orth ls_start = P1 + i * ls_delta;
         clipper::Coord_orth ls_end   = ls_start + ls_delta;
         Cartesian ls_start_c(ls_start);
         Cartesian ls_end_c(ls_end);

         loop_line_segments.push_back(CartesianPair(ls_start_c, ls_end_c));
      }
   }
   return std::pair<bool, std::vector<coot::CartesianPair>>(needs_CA_CA_badness_spots, loop_line_segments);
}

} // namespace coot