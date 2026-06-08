#pragma once

#include <string>
#include <gemmi/model.hpp>
#include "geometry/protein-geometry.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"
#include "clipper/core/coords.h"
#include <cmath>

namespace coot {
namespace restraints_container_t {

static bool dictionary_name_matches_coords_resname(const std::string &comp_id,
                                                    const std::string &resname) {
   std::string r = resname;
   if (r.length() > 2)
      if (r[2] == ' ')
         r = resname.substr(0, 2);
   return (r == comp_id);
}

inline std::string find_peptide_link_type_ng_gemmi(
    const gemmi::Residue& res_1,
    const gemmi::Residue& res_2,
    const coot::protein_geometry &geom) {

   std::string link_type;

   std::string t1;
   std::string t2;

   std::string residue_type_1 = res_1.name;
   std::string residue_type_2 = res_2.name;

   for (unsigned int idr = 0; idr < geom.size(); idr++) {
      if (dictionary_name_matches_coords_resname(geom.three_letter_code(idr), residue_type_1)) {
         t1 = geom[idr].second.residue_info.group;
         break;
      }
   }
   for (unsigned int idr = 0; idr < geom.size(); idr++) {
      if (dictionary_name_matches_coords_resname(geom.three_letter_code(idr), residue_type_2)) {
         t2 = geom[idr].second.residue_info.group;
         break;
      }
   }

   if (t1 == "L-peptide" || t1 == "D-peptide" || t1 == "M-peptide" || t1 == "P-peptide" || t1 == "peptide") {
      if (t2 == "L-peptide" || t2 == "D-peptide" || t2 == "M-peptide" || t2 == "P-peptide" || t2 == "peptide") {
         if (residue_type_2 == "PRO" || residue_type_2 == "HYP") {
            link_type = "PTRANS";
         } else {
            link_type = "TRANS";
         }
      }
   }

   const gemmi::Atom* ca_r_1 = nullptr;
   const gemmi::Atom* c_r_1  = nullptr;
   const gemmi::Atom* n_r_2  = nullptr;
   const gemmi::Atom* ca_r_2 = nullptr;

   for (const auto& at : res_1.atoms) {
      std::string at_name(at.name);
      if (at_name == "C")
         c_r_1 = &at;
      if (at_name == "CA")
         ca_r_1 = &at;
   }
   for (const auto& at : res_2.atoms) {
      std::string at_name(at.name);
      if (at_name == "N")
         n_r_2 = &at;
      if (at_name == "CA")
         ca_r_2 = &at;
   }

   if (ca_r_1 && c_r_1 && n_r_2 && ca_r_2) {
      clipper::Coord_orth pt_1 = coot::co_gemmi(ca_r_1);
      clipper::Coord_orth pt_2 = coot::co_gemmi(c_r_1);
      clipper::Coord_orth pt_3 = coot::co_gemmi(n_r_2);
      clipper::Coord_orth pt_4 = coot::co_gemmi(ca_r_2);
      double torsion = clipper::Coord_orth::torsion(pt_1, pt_2, pt_3, pt_4);
      if (torsion > M_PI)
         torsion -= 2.0 * M_PI;
      double ninety_degrees = 0.5 * M_PI;
      if (torsion < ninety_degrees && torsion > -ninety_degrees) {
         if (link_type == "TRANS") link_type = "CIS";
         if (link_type == "PTRANS") link_type = "PCIS";
      }
   }

   return link_type;
}

} // namespace restraints_container_t
} // namespace coot
