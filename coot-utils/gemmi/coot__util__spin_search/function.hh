#pragma once

#include <utility>
#include <vector>
#include <iostream>
#include <string>

#include <clipper/core/xmap.h>
#include "clipper/core/coords.h"
#include "clipper/core/clipper_util.h"

#include <gemmi/model.hpp>

#include "geometry/residue-and-atom-specs.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/coot-map-utils.hh"

namespace coot {
namespace util {

std::pair<float, float> spin_search_gemmi(
    const clipper::Xmap<float>& xmap,
    gemmi::Chain* chain,
    gemmi::Residue* residue,
    coot::torsion tors
) {

   float best_ori = -1111.1f;
   float torsion_relative_to_N = -1111.1f;

   // Inline matching_atoms for gemmi types
   // Match each of the 4 atom specs against atoms in the residue
   std::vector<gemmi::Atom*> match_atoms;
   
   struct spec_ref {
       const coot::atom_spec_t* spec;
   };
   spec_ref refs[] = {
       {&tors.atom_1.second},
       {&tors.atom_2.second},
       {&tors.atom_3.second},
       {&tors.atom_4.second}
   };

   for (int si = 0; si < 4; si++) {
       const coot::atom_spec_t* spec = refs[si].spec;
       gemmi::Atom* found = nullptr;
       for (const gemmi::Atom& atom : residue->atoms) {
           // Match atom name: strip both to raw names (gemmi atom.name is unpadded,
           // spec->atom_name is 4-char padded from PDB columns)
           std::string gemmi_name = atom.name;
           std::string spec_name = spec->atom_name;
           // Remove trailing spaces from spec name
           while (!spec_name.empty() && spec_name.back() == ' ') spec_name.pop_back();
           // Remove leading spaces from spec name  
           auto first_non_space = spec_name.find_first_not_of(' ');
           if (first_non_space != std::string::npos) spec_name = spec_name.substr(first_non_space);
           else spec_name.clear();
           if (spec_name != gemmi_name) continue;
           // Match alt_conf (empty string = no alt, gemmi uses space char)
           if (!spec->alt_conf.empty() && spec->alt_conf[0] != atom.altloc) continue;
           // Match residue number
           if (spec->res_no != residue->seqid.num.value) continue;
           // Match insertion code (empty string = no ins code, gemmi uses space)
           if (!spec->ins_code.empty() && spec->ins_code[0] != residue->seqid.icode) continue;
           // Match chain
           if (chain && !spec->chain_id.empty() && spec->chain_id != chain->name) continue;
           found = const_cast<gemmi::Atom*>(&atom);
           break;
       }
       if (found) {
           match_atoms.push_back(found);
       }
   }

   if (match_atoms.size() != 4) {
      std::cout << "ERROR:: not all atoms for torsion found in residue!" << std::endl;
      std::cout << "        (found " << match_atoms.size() << " atoms.)" << std::endl;
   } else {

      clipper::Coord_orth pa1(match_atoms[0]->pos.x, match_atoms[0]->pos.y, match_atoms[0]->pos.z);
      clipper::Coord_orth pa2(match_atoms[1]->pos.x, match_atoms[1]->pos.y, match_atoms[1]->pos.z);
      clipper::Coord_orth pa3(match_atoms[2]->pos.x, match_atoms[2]->pos.y, match_atoms[2]->pos.z);
      clipper::Coord_orth pa4(match_atoms[3]->pos.x, match_atoms[3]->pos.y, match_atoms[3]->pos.z);

      float best_d = -99999999.9f;
      clipper::Coord_orth best_pos;
      for (double theta = 0; theta <= 360; theta += 3.0) {

         clipper::Coord_orth dir   = pa3 - pa2;
         clipper::Coord_orth pos   = pa4;
         clipper::Coord_orth shift = pa3;
         clipper::Coord_orth co = coot::util::rotate_around_vector(dir, pos, shift, clipper::Util::d2rad(theta));
         float this_d = coot::util::density_at_point(xmap, co);
         if (this_d > best_d) {
            best_d = this_d;
            best_ori = static_cast<float>(theta);
            best_pos = co;
         }
      }
      double ct = clipper::Coord_orth::torsion(pa1, pa2, pa3, best_pos);
      torsion_relative_to_N = static_cast<float>(clipper::Util::rad2d(ct));
   }

   std::pair<float, float> p(best_ori, torsion_relative_to_N);
   return p;
}

} // namespace util
} // namespace coot
