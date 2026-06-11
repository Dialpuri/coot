#pragma once

#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include <clipper/core/clipper_util.h>
#include <map>
#include <string>
#include <cmath>

namespace coot {

   // gemmi version of atom_quad using gemmi::Atom* instead of mmdb::Atom*
   struct atom_quad_gemmi {
      gemmi::Atom *atom_1; // N
      gemmi::Atom *atom_2; // CA
      gemmi::Atom *atom_3; // C
      gemmi::Atom *atom_4; // CB

      atom_quad_gemmi() : atom_1(nullptr), atom_2(nullptr), atom_3(nullptr), atom_4(nullptr) {}

      bool filled_p() const {
         return atom_1 != nullptr && atom_2 != nullptr && atom_3 != nullptr && atom_4 != nullptr;
      }
   };

   // gemmi-compatible c_beta_deviation using gemmi::Atom*
   struct c_beta_deviation_gemmi {
      gemmi::Atom *at;
      clipper::Coord_orth pos_ideal;
      double dist;

      c_beta_deviation_gemmi(gemmi::Atom *a, const clipper::Coord_orth &p, double d)
         : at(a), pos_ideal(p), dist(d) {}

      c_beta_deviation_gemmi()
         : at(nullptr), pos_ideal(clipper::Coord_orth(-1,-1,-1)), dist(0.0) {}
   };

   inline clipper::Coord_orth make_CB_ideal_pos_gemmi(const atom_quad_gemmi &q, const std::string &res_name) {
      const double l = 1.52;
      double a1 = 0.0, a2 = 0.0;
      double t1 = 120.0, t2 = -120.0;

      if (res_name == "ALA") a1 = clipper::Util::d2rad(109.344);
      if (res_name == "CYS") a1 = clipper::Util::d2rad(109.344);
      if (res_name == "ASP") a1 = clipper::Util::d2rad(109.344);
      if (res_name == "GLU") a1 = clipper::Util::d2rad(109.344);
      if (res_name == "PHE") a1 = clipper::Util::d2rad(109.344);
      if (res_name == "GLY") a1 = clipper::Util::d2rad(109.344);
      if (res_name == "HIS") a1 = clipper::Util::d2rad(109.344);
      if (res_name == "ILE") a1 = clipper::Util::d2rad(109.344);
      if (res_name == "LYS") a1 = clipper::Util::d2rad(109.344);
      if (res_name == "LEU") a1 = clipper::Util::d2rad(109.344);
      if (res_name == "MET") a1 = clipper::Util::d2rad(109.344);
      if (res_name == "ASN") a1 = clipper::Util::d2rad(109.344);
      if (res_name == "PRO") a1 = clipper::Util::d2rad(109.344);
      if (res_name == "GLN") a1 = clipper::Util::d2rad(109.344);
      if (res_name == "ARG") a1 = clipper::Util::d2rad(109.344);
      if (res_name == "SER") a1 = clipper::Util::d2rad(109.344);
      if (res_name == "THR") a1 = clipper::Util::d2rad(109.344);
      if (res_name == "VAL") a1 = clipper::Util::d2rad(109.344);
      if (res_name == "TRP") a1 = clipper::Util::d2rad(109.344);
      if (res_name == "TYR") a1 = clipper::Util::d2rad(109.344);

      a2 = a1;

      if (res_name == "CYS") a2 = clipper::Util::d2rad(109.344);
      if (res_name == "ASP") a2 = clipper::Util::d2rad(109.344);
      if (res_name == "GLU") a2 = clipper::Util::d2rad(109.344);
      if (res_name == "PHE") a2 = clipper::Util::d2rad(109.344);
      if (res_name == "HIS") a2 = clipper::Util::d2rad(109.344);
      if (res_name == "ILE") a2 = clipper::Util::d2rad(109.344);
      if (res_name == "LYS") a2 = clipper::Util::d2rad(109.344);
      if (res_name == "LEU") a2 = clipper::Util::d2rad(109.344);
      if (res_name == "MET") a2 = clipper::Util::d2rad(109.344);
      if (res_name == "PRO") a2 = clipper::Util::d2rad(110.031);
      if (res_name == "GLN") a2 = clipper::Util::d2rad(111.037);
      if (res_name == "ARG") a2 = clipper::Util::d2rad(111.037);
      if (res_name == "SER") a2 = clipper::Util::d2rad(111.379);
      if (res_name == "THR") a2 = clipper::Util::d2rad(111.511);
      if (res_name == "VAL") a2 = clipper::Util::d2rad(111.388);
      if (res_name == "TRP") a2 = clipper::Util::d2rad(111.644);
      if (res_name == "TYR") a2 = clipper::Util::d2rad(111.331);
      if (res_name == "PRO") {
         t1 = clipper::Util::d2rad(115.20);
         t2 = clipper::Util::d2rad(-119.80);
      }

      clipper::Coord_orth pt_1(q.atom_1->pos.x, q.atom_1->pos.y, q.atom_1->pos.z);
      clipper::Coord_orth pt_2(q.atom_2->pos.x, q.atom_2->pos.y, q.atom_2->pos.z);
      clipper::Coord_orth pt_3(q.atom_3->pos.x, q.atom_3->pos.y, q.atom_3->pos.z);

      clipper::Coord_orth pt_trial1 = clipper::Coord_orth(pt_1, pt_3, pt_2, l, a2, t1);
      clipper::Coord_orth pt_trial2 = clipper::Coord_orth(pt_3, pt_1, pt_2, l, a1, t2);

      return 0.5 * (pt_trial1 + pt_trial2);
   }

   inline std::string trim_atom_name(const std::string &s) {
      size_t start = s.find_first_not_of(" \t\r\n");
      if (start == std::string::npos) return "";
      size_t end = s.find_last_not_of(" \t\r\n");
      return s.substr(start, end - start + 1);
   }

   inline
   std::map<std::string, c_beta_deviation_gemmi>
   get_c_beta_deviations_gemmi(gemmi::Residue *residue_p) {

      std::map<std::string, c_beta_deviation_gemmi> m;
      std::string res_name(residue_p->name);

      std::map<std::string, atom_quad_gemmi> alt_conf_map;

      for (gemmi::Atom &at : residue_p->atoms) {
         std::string atom_name = trim_atom_name(at.name);
         std::string alt_conf;
         if (at.altloc != ' ' && at.altloc != '\0') {
            alt_conf = std::string(1, at.altloc);
         }
         // else alt_conf stays "" (empty string), matching MMDB behavior
         if (atom_name == "N")  alt_conf_map[alt_conf].atom_1 = &at;
         if (atom_name == "CA") alt_conf_map[alt_conf].atom_2 = &at;
         if (atom_name == "C")  alt_conf_map[alt_conf].atom_3 = &at;
         if (atom_name == "CB") alt_conf_map[alt_conf].atom_4 = &at;
      }

      for (auto const &kv : alt_conf_map) {
         const atom_quad_gemmi &q = kv.second;
         if (q.filled_p()) {
            clipper::Coord_orth CB_real_pos(q.atom_4->pos.x, q.atom_4->pos.y, q.atom_4->pos.z);
            clipper::Coord_orth CB_ideal_pos = make_CB_ideal_pos_gemmi(q, res_name);
            double dsqrd = (CB_ideal_pos - CB_real_pos).lengthsq();
            double d = sqrt(dsqrd);
            m[kv.first] = c_beta_deviation_gemmi(q.atom_4, CB_ideal_pos, d);
         }
      }

      return m;
   }

} // namespace coot
