#include "function.hh"
#include <sstream>
#include <clipper/core/coords.h>

namespace coot {
namespace dipole {

void init_gemmi(
    dipole_state_t& result,
    std::vector<std::pair<coot::dictionary_residue_restraints_t, gemmi::CRA>> dict_res_pairs
) {
   result.dipole_is_good_flag = false;

   // Find the centre of the residue
   double sum_x = 0;
   double sum_y = 0;
   double sum_z = 0;
   int n_points = 0;

   for (unsigned int ires = 0; ires < dict_res_pairs.size(); ires++) {
      gemmi::CRA cra = dict_res_pairs[ires].second;
      coot::dictionary_residue_restraints_t rest = dict_res_pairs[ires].first;

      if (cra.residue == nullptr) continue;

      for (gemmi::Atom& atom : cra.residue->atoms) {
         sum_x += atom.pos.x;
         sum_y += atom.pos.y;
         sum_z += atom.pos.z;
         n_points++;
      }
   }

   if (n_points == 0) {
      std::string mess = "No atoms in ";
      mess += std::to_string(dict_res_pairs.size());
      mess += " residue";
      if (dict_res_pairs.size() != 1)
         mess += "s";
      mess += " ";
      for (unsigned int i = 0; i < dict_res_pairs.size(); i++) {
         gemmi::CRA cra = dict_res_pairs[i].second;
         if (cra.chain != nullptr)
            mess += cra.chain->name;
         mess += " ";
         if (cra.residue != nullptr)
            mess += std::to_string(cra.residue->seqid.num.value);
         mess += ", ";
      }
      throw std::runtime_error(mess);
   }

   double multiplier = 1.0 / double(n_points);

   // class member:
   result.residue_centre = clipper::Coord_orth(sum_x * multiplier,
                                               sum_y * multiplier,
                                               sum_z * multiplier);

   // Convert dict_res_pairs for charged_atoms_gemmi call
   // charged_atoms_gemmi expects: std::vector<std::pair<coot::dictionary_residue_restraints_t, gemmi::Residue*>>
   std::vector<std::pair<coot::dictionary_residue_restraints_t, gemmi::Residue*>> charged_atoms_input;
   for (auto& p : dict_res_pairs) {
      charged_atoms_input.push_back(std::make_pair(p.first, p.second.residue));
   }

   std::vector<std::pair<gemmi::Atom*, float>> charged_ats =
       coot::dipole::charged_atoms_gemmi(charged_atoms_input);

   std::vector<std::pair<float, clipper::Coord_orth>> charged_points(charged_ats.size());
   for (unsigned int i = 0; i < charged_ats.size(); i++) {
      clipper::Coord_orth p(charged_ats[i].first->pos.x,
                            charged_ats[i].first->pos.y,
                            charged_ats[i].first->pos.z);
      charged_points[i] =
         std::pair<float, clipper::Coord_orth>(charged_ats[i].second, p);
   }

   clipper::Coord_orth dip_val(0, 0, 0);
   for (unsigned int ii = 0; ii < charged_points.size(); ii++) {
      clipper::Coord_orth scaled(charged_points[ii].second - result.residue_centre);
      scaled = charged_points[ii].first * scaled;
      dip_val += scaled;
      result.dipole_is_good_flag = true;
   }

   if (!result.dipole_is_good_flag) {
      std::string mess = "Dipole is not good for ";
      mess += std::to_string(dict_res_pairs.size());
      mess += " residue";
      if (dict_res_pairs.size() != 1)
         mess += "s";
      mess += " ";
      for (unsigned int i = 0; i < dict_res_pairs.size(); i++) {
         gemmi::CRA cra = dict_res_pairs[i].second;
         if (cra.chain != nullptr)
            mess += cra.chain->name;
         mess += " ";
         if (cra.residue != nullptr) {
            mess += std::to_string(cra.residue->seqid.num.value);
            mess += " ";
            mess += cra.residue->name;
         }
         mess += ", ";
      }
      throw std::runtime_error(mess);
   }

   result.dipole_ = dip_val;
}

} // namespace dipole
} // namespace coot