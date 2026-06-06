#pragma once

#include <vector>
#include <string>
#include <iostream>

#include <gemmi/model.hpp>

#include "mini-mol/mini-mol.hh"
#include "clipper/core/clipper_types.h"

namespace coot {
namespace high_res {

inline coot::minimol::molecule
filter_on_groups_gemmi(const std::vector<std::vector<int>>& groups,
                       std::vector<gemmi::CRA> atom_selection) {

   coot::minimol::molecule m;
   for (unsigned int igroup = 0; igroup < groups.size(); igroup++) {

      std::cout << "group " << igroup << " has " << groups[igroup].size()
                << " members" << std::endl;

      // Find the average position in the group
      clipper::Coord_orth sum(0.0, 0.0, 0.0);
      int n_grp_ats = groups[igroup].size();
      for (int iat = 0; iat < n_grp_ats; iat++) {
         gemmi::CRA cra = atom_selection[groups[igroup][iat]];
         clipper::Coord_orth pt(cra.atom->pos.x, cra.atom->pos.y, cra.atom->pos.z);
         sum += pt;
      }
      double frac = 1.0 / double(n_grp_ats);
      clipper::Coord_orth new_sum(0.0, 0.0, 0.0);
      clipper::Coord_orth mean_pt(sum.x()*frac, sum.y()*frac, sum.z()*frac);
      int n_group_within_lim = 0;
      for (int iat = 0; iat < n_grp_ats; iat++) {
         gemmi::CRA cra = atom_selection[groups[igroup][iat]];
         clipper::Coord_orth pt(cra.atom->pos.x, cra.atom->pos.y, cra.atom->pos.z);
         n_group_within_lim++;
         new_sum += pt;
      }
      if (groups[igroup].size() > 10) {
         for (unsigned int iat = 0; iat < groups[igroup].size(); iat++) {
            gemmi::CRA cra = atom_selection[groups[igroup][iat]];
            std::cout << "   " << cra.atom->name << std::endl;
            std::cout << "   (set-rotation-centre "
                      << cra.atom->pos.x << " "
                      << cra.atom->pos.y << " "
                      << cra.atom->pos.z << ")" << std::endl;
         }
      }

      if (n_group_within_lim > 0) {
         frac = 1.0 / double(n_group_within_lim);
         clipper::Coord_orth av_pt(new_sum.x()*frac,
                                   new_sum.y()*frac,
                                   new_sum.z()*frac);
         gemmi::CRA spec_cra = atom_selection[groups[igroup][0]];
         std::string atom_name(spec_cra.atom->name);
         std::string atom_element(spec_cra.atom->element.name());
         int resno = spec_cra.residue->seqid.num.value;
         std::string chain_id(spec_cra.chain->name);
         int ifrag = m.fragment_for_chain(chain_id);
         coot::minimol::atom atom(atom_name, atom_element, av_pt, "", 30.0);
         m[ifrag][resno].name = "ALA";
         m[ifrag][resno].seqnum = resno;
         m[ifrag][resno].addatom(atom);
      } else {
         gemmi::CRA cra = atom_selection[groups[igroup][0]];
         std::cout << "OOps! No unfiltered atoms left for "
                   << cra.atom->name << std::endl;
      }
   }
   return m;
}

} // high_res
} // coot