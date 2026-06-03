#pragma once

#include <map>
#include <string>
#include <cmath>
#include <gemmi/model.hpp>
#include <clipper/core/coords.h>

namespace coot {

   class c_beta_deviation_t {
   public:
      c_beta_deviation_t(gemmi::Atom *at_in, const clipper::Coord_orth &pos_ideal_in, double dist_in) :
         at(at_in), pos_ideal(pos_ideal_in), dist(dist_in) { }
      c_beta_deviation_t() : at(nullptr), pos_ideal(clipper::Coord_orth(-1,-1,-1)), dist(0.0) {}
      gemmi::Atom *at;
      clipper::Coord_orth pos_ideal;
      double dist;
   };

   // Local replacement for coot::atom_quad — original uses mmdb::Atom*
   struct atom_quad_gemmi {
      gemmi::Atom *atom_1;
      gemmi::Atom *atom_2;
      gemmi::Atom *atom_3;
      gemmi::Atom *atom_4;
      atom_quad_gemmi() : atom_1(nullptr), atom_2(nullptr), atom_3(nullptr), atom_4(nullptr) {}
      bool filled_p() const {
         return atom_1 != nullptr && atom_2 != nullptr && atom_3 != nullptr && atom_4 != nullptr;
      }
   };

   inline
   clipper::Coord_orth
   make_CB_ideal_pos_gemmi(const atom_quad_gemmi &q, const std::string &res_name)
   {
      bool is_PRO = false;
      if (res_name == std::string("PRO"))
         is_PRO = true;

      clipper::Coord_orth pt_1(q.atom_1->pos.x, q.atom_1->pos.y, q.atom_1->pos.z);
      clipper::Coord_orth pt_2(q.atom_2->pos.x, q.atom_2->pos.y, q.atom_2->pos.z);
      clipper::Coord_orth pt_3(q.atom_3->pos.x, q.atom_3->pos.y, q.atom_3->pos.z);

      double l = 1.53;
      double a1 = clipper::Util::d2rad(111.0);
      double a2 = clipper::Util::d2rad(111.7);
      double t1 = clipper::Util::d2rad(122.9);
      double t2 = clipper::Util::d2rad(-122.6);

      if (res_name == "ALA") l = 1.509;
      if (res_name == "ASP") l = 1.531;
      if (res_name == "ASN") l = 1.531;
      if (res_name == "CYS") l = 1.524;
      if (res_name == "GLU") l = 1.530;
      if (res_name == "PHE") l = 1.531;
      if (res_name == "HIS") l = 1.535;
      if (res_name == "ILE") l = 1.542;
      if (res_name == "LYS") l = 1.532;
      if (res_name == "LEU") l = 1.532;
      if (res_name == "MET") l = 1.532;
      if (res_name == "PRO") l = 1.534;
      if (res_name == "GLN") l = 1.530;
      if (res_name == "ARG") l = 1.532;
      if (res_name == "SER") l = 1.507;
      if (res_name == "THR") l = 1.534;
      if (res_name == "VAL") l = 1.541;
      if (res_name == "TRP") l = 1.534;
      if (res_name == "TYR") l = 1.531;

      if (res_name == "ALA") a1 = clipper::Util::d2rad(109.912);
      if (res_name == "ASP") a1 = clipper::Util::d2rad(111.338);
      if (res_name == "ASN") a1 = clipper::Util::d2rad(111.766);
      if (res_name == "CYS") a1 = clipper::Util::d2rad(110.827);
      if (res_name == "GLU") a1 = clipper::Util::d2rad(110.374);
      if (res_name == "PHE") a1 = clipper::Util::d2rad(110.494);
      if (res_name == "HIS") a1 = clipper::Util::d2rad(110.437);
      if (res_name == "ILE") a1 = clipper::Util::d2rad(110.820);
      if (res_name == "LYS") a1 = clipper::Util::d2rad(110.374);
      if (res_name == "LEU") a1 = clipper::Util::d2rad(108.955);
      if (res_name == "MET") a1 = clipper::Util::d2rad(110.906);
      if (res_name == "PRO") a1 = clipper::Util::d2rad(103.430);
      if (res_name == "GLN") a1 = clipper::Util::d2rad(110.374);
      if (res_name == "ARG") a1 = clipper::Util::d2rad(110.374);
      if (res_name == "SER") a1 = clipper::Util::d2rad(110.990);
      if (res_name == "THR") a1 = clipper::Util::d2rad(111.125);
      if (res_name == "VAL") a1 = clipper::Util::d2rad(111.441);
      if (res_name == "TRP") a1 = clipper::Util::d2rad(110.562);
      if (res_name == "TYR") a1 = clipper::Util::d2rad(110.494);

      if (res_name == "ALA") a2 = clipper::Util::d2rad(111.490);
      if (res_name == "ASP") a2 = clipper::Util::d2rad(111.804);
      if (res_name == "ASN") a2 = clipper::Util::d2rad(111.540);
      if (res_name == "CYS") a2 = clipper::Util::d2rad(109.612);
      if (res_name == "GLU") a2 = clipper::Util::d2rad(111.037);
      if (res_name == "PHE") a2 = clipper::Util::d2rad(111.331);
      if (res_name == "HIS") a2 = clipper::Util::d2rad(112.128);
      if (res_name == "ILE") a2 = clipper::Util::d2rad(111.764);
      if (res_name == "LYS") a2 = clipper::Util::d2rad(111.037);
      if (res_name == "LEU") a2 = clipper::Util::d2rad(111.075);
      if (res_name == "MET") a2 = clipper::Util::d2rad(109.344);
      if (res_name == "PRO") a2 = clipper::Util::d2rad(110.031);
      if (res_name == "GLN") a2 = clipper::Util::d2rad(111.037);
      if (res_name == "ARG") a2 = clipper::Util::d2rad(111.037);
      if (res_name == "SER") a2 = clipper::Util::d2rad(111.379);
      if (res_name == "THR") a2 = clipper::Util::d2rad(111.511);
      if (res_name == "VAL") a2 = clipper::Util::d2rad(111.388);
      if (res_name == "TRP") a2 = clipper::Util::d2rad(111.644);
      if (res_name == "TYR") a2 = clipper::Util::d2rad(111.331);
      if (is_PRO) {
         t1 = clipper::Util::d2rad(115.20);
         t2 = clipper::Util::d2rad(-119.80);
      }
      clipper::Coord_orth pt_trial1 = clipper::Coord_orth(pt_1, pt_3, pt_2, l, a2, t1);
      clipper::Coord_orth pt_trial2 = clipper::Coord_orth(pt_3, pt_1, pt_2, l, a1, t2);

      clipper::Coord_orth pt_trial = 0.5 * (pt_trial1 + pt_trial2);

      return pt_trial;
   }

   // Trim leading and trailing whitespace from a string
   inline std::string trim(const std::string &s) {
      size_t start = s.find_first_not_of(" \t\r\n");
      if (start == std::string::npos) return "";
      size_t end = s.find_last_not_of(" \t\r\n");
      return s.substr(start, end - start + 1);
   }

   inline
   std::map<std::string, coot::c_beta_deviation_t>
   get_c_beta_deviations_residue_gemmi(gemmi::Residue &residue)
   {
      std::map<std::string, c_beta_deviation_t> m;
      std::string res_name(residue.name);

      std::map<std::string, atom_quad_gemmi> alt_conf_map;

      for (gemmi::Atom &at : residue.atoms) {
         std::string atom_name = trim(at.name);
         std::string alt_conf(1, at.altloc);
         if (atom_name == "N") alt_conf_map[alt_conf].atom_1 = &at;
         if (atom_name == "CA") alt_conf_map[alt_conf].atom_2 = &at;
         if (atom_name == "C") alt_conf_map[alt_conf].atom_3 = &at;
         if (atom_name == "CB") alt_conf_map[alt_conf].atom_4 = &at;
      }

      for (auto const &kv : alt_conf_map) {
         const atom_quad_gemmi &q = kv.second;
         if (q.filled_p()) {
            clipper::Coord_orth CB_real_pos(q.atom_4->pos.x, q.atom_4->pos.y, q.atom_4->pos.z);
            clipper::Coord_orth CB_ideal_pos = make_CB_ideal_pos_gemmi(q, res_name);
            double dsqrd = (CB_ideal_pos - CB_real_pos).lengthsq();
            double d = sqrt(dsqrd);
            c_beta_deviation_t cbd(q.atom_4, CB_ideal_pos, d);
            m[kv.first] = cbd;
         }
      }

      return m;
   }

   inline
   std::map<const gemmi::Residue*, std::map<std::string, c_beta_deviation_t>>
   get_c_beta_deviations_gemmi(gemmi::Structure &st)
   {
      std::map<const gemmi::Residue*, std::map<std::string, c_beta_deviation_t>> m;

      if (!st.models.empty()) {
         for (gemmi::Chain &chain : st.models[0].chains) {
            for (gemmi::Residue &residue : chain.residues) {
               std::map<std::string, c_beta_deviation_t> cbdm = get_c_beta_deviations_residue_gemmi(residue);
               if (!cbdm.empty()) {
                  m[&residue] = cbdm;
               }
            }
         }
      }
      return m;
   }

} // namespace coot