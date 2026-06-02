#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <string>
#include <cmath>
#include <clipper/clipper.h>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__calc_cablam/gemmi/function.hh"

namespace coot {

class cablam_like_geometry_stats_t_gemmi {
public:
   gemmi::Residue *residue;
   std::string chain_id;
   double dp_prev_to_mid;
   double dp_next_to_mid;
   double dist_proj_point_prev_to_next;

   cablam_like_geometry_stats_t_gemmi(const cablam_markup_t_gemmi &cm, const std::string& chain_id);
};

inline cablam_like_geometry_stats_t_gemmi::cablam_like_geometry_stats_t_gemmi(
    const cablam_markup_t_gemmi &cm, const std::string& chain_id) {

   auto v_prev = cm.O_prev_pos - cm.CA_proj_point_prev;
   auto v_this = cm.O_this_pos - cm.CA_proj_point_this;
   auto v_next = cm.O_next_pos - cm.CA_proj_point_next;

   double dd = (cm.CA_proj_point_next - cm.CA_proj_point_prev).lengthsq();
   double d = std::sqrt(dd);
   double a1 = clipper::Coord_orth::dot(v_prev, v_this);
   double a2 = clipper::Coord_orth::dot(v_next, v_this);

   residue = cm.residue;
   this->chain_id = chain_id;
   dp_prev_to_mid = a1;
   dp_next_to_mid = a2;
   dist_proj_point_prev_to_next = d;
}

std::vector<cablam_like_geometry_stats_t_gemmi>
get_cablam_like_geometry_stats_gemmi(gemmi::Structure& mol) {

   std::vector<cablam_like_geometry_stats_t_gemmi> v;

   if (mol.models.empty()) return v;

   gemmi::Model& model = mol.models[0];

   for (gemmi::Chain& chain : model.chains) {
      int n_res = static_cast<int>(chain.residues.size());
      int n_res_max = n_res - 2;
      for (int ires = 1; ires < n_res_max; ires++) {
         gemmi::Residue* residue_p = &chain.residues[ires];
         auto cm = calc_cablam_gemmi(&chain, residue_p, ires, 0.0);
         if (cm.residue) {
            v.push_back(cablam_like_geometry_stats_t_gemmi(cm, chain.name));
         }
      }
   }

   return v;
}

} // namespace coot