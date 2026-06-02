#pragma once

#include <vector>
#include <gemmi/model.hpp>
#include <clipper/core/coords.h>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__compare_to_helix__188839/gemmi/function.hh"

namespace coot {

std::vector<gemmi::Residue *>
like_a_helix_gemmi(const std::vector<gemmi::Residue *> &selected_residues) {

   std::vector<gemmi::Residue *> helical_residues;

   int n_sel = static_cast<int>(selected_residues.size());

   // ideal helix atom positions (N CA C O backbone for each residue in a turn)
   std::vector<clipper::Coord_orth> ref_pos;
   ref_pos.push_back(clipper::Coord_orth(0.370,   3.598,   2.351));
   ref_pos.push_back(clipper::Coord_orth(1.039,   3.821,   3.641));
   ref_pos.push_back(clipper::Coord_orth(2.491,   3.342,   3.566));
   ref_pos.push_back(clipper::Coord_orth(3.418,   4.049,   3.988));
   ref_pos.push_back(clipper::Coord_orth(2.634,   2.146,   3.027));
   ref_pos.push_back(clipper::Coord_orth(3.942,   1.496,   2.858));
   ref_pos.push_back(clipper::Coord_orth(4.862,   2.382,   2.015));
   ref_pos.push_back(clipper::Coord_orth(6.030,   2.606,   2.366));
   ref_pos.push_back(clipper::Coord_orth(4.295,   2.858,   0.922));
   ref_pos.push_back(clipper::Coord_orth(4.999,   3.730,  -0.029));
   ref_pos.push_back(clipper::Coord_orth(5.506,   4.984,   0.686));
   ref_pos.push_back(clipper::Coord_orth(6.667,   5.386,   0.524));
   ref_pos.push_back(clipper::Coord_orth(4.606,   5.561,   1.460));
   ref_pos.push_back(clipper::Coord_orth(4.881,   6.777,   2.238));
   ref_pos.push_back(clipper::Coord_orth(6.063,   6.539,   3.181));
   ref_pos.push_back(clipper::Coord_orth(6.987,   7.362,   3.267));

   if (n_sel > 3) {

      for (int istart = 0; istart < (n_sel - 4); istart++) {
         std::vector<gemmi::Residue *> test_residues;
         if ((istart + 4) < n_sel) {
            for (int i_4 = istart; i_4 < (istart + 4); i_4++)
               test_residues.push_back(selected_residues[i_4]);

            if (test_residues.size() == 4) {
               helical_results_t helicals =
                   compare_to_helix_gemmi(test_residues, ref_pos);
               if (helicals.is_alpha_helix_like)
                  helical_residues.push_back(test_residues[1]);
            }
         }
      }
   }
   return helical_residues;
}

} // namespace coot