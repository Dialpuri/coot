/* coot-utils/coot-coord-utils-nucleotides.cc
 * 
 * Copyright 2008, 2009 by The University of Oxford
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#include <optional>
#include <stdexcept>

#include "clipper/core/coords.h"
#include "geometry/residue-and-atom-specs.hh"
#include "utils/coot-utils.hh"
#include "coot-coord-utils.hh"

 
// Throw an exception if it is not possible to generate pucker info
// 
coot::pucker_analysis_info_t::pucker_analysis_info_t(mmdb::Residue *res_p,
                                                     std::string altconf_in) {

   auto get_base_lsq_plane = [] (const std::vector<clipper::Coord_orth> &coords) -> 
      std::optional<lsq_plane_info_t> {

      if (coords.size() < 3) return std::nullopt;
      lsq_plane_info_t lsq_plane(coords);
      return lsq_plane;
   };

   out_of_plane_distance = 0.0;
   plane_distortion = 0.0;
   C1_prime = nullptr;
   N1_or_9  = nullptr;
   // The atoms are in the following order C1' C2' C3' C4' O4
   //
   altconf = altconf_in; // save for phosphate distance (if needed).

   assign_base_atom_coords(res_p); // and C1_prime and N1_or_9 if possible
   std::optional<lsq_plane_info_t> lsq_plane = get_base_lsq_plane(base_atoms_coords);

   if (lsq_plane.has_value()) {

      // store the geometry
      markup_info.base_ring_centre = lsq_plane.value().centre();
      markup_info.base_ring_normal = lsq_plane.value().normal();

      std::vector<mmdb::Atom *> ribose_atoms(5, nullptr); // ribose atoms
      std::vector<coot::pucker_analysis_info_t::PUCKERED_ATOM_T> possible_puckers;
      possible_puckers.push_back(coot::pucker_analysis_info_t::C1_PRIME);
      possible_puckers.push_back(coot::pucker_analysis_info_t::C2_PRIME);
      possible_puckers.push_back(coot::pucker_analysis_info_t::C3_PRIME);
      possible_puckers.push_back(coot::pucker_analysis_info_t::C4_PRIME);
      possible_puckers.push_back(coot::pucker_analysis_info_t::O4_PRIME);

      mmdb::PPAtom residue_atoms = NULL;
      int n_residue_atoms = 0;
      res_p->GetAtomTable(residue_atoms, n_residue_atoms);
      // find the phosphorus atom
      for (int i=0; i<n_residue_atoms; i++) {
         mmdb::Atom *atm = residue_atoms[i];
         if (! atm->isTer()) {
            std::string atm_name(atm->name);
            std::string alt_name(atm->altLoc);
            if (altconf == alt_name) {
               if (atm_name == " P  ") { // PDBv3 FIXME
                  clipper::Coord_orth p(atm->x, atm->y, atm->z);
                  markup_info.phosphorus_position = p;
                  markup_info.projected_point = lsq_plane.value().projected_point(p);
               }
            }
         }
      }
      // find the ribose atoms
      res_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int i=0; i<n_residue_atoms; i++) {
         std::string atm_name(residue_atoms[i]->name);
         std::string alt_name(residue_atoms[i]->altLoc);
         if (altconf == alt_name) {
            if (atm_name == " C1*") ribose_atoms[0] = residue_atoms[i];
            if (atm_name == " C1'") ribose_atoms[0] = residue_atoms[i];
            if (atm_name == " C2*") ribose_atoms[1] = residue_atoms[i];
            if (atm_name == " C2'") ribose_atoms[1] = residue_atoms[i];
            if (atm_name == " C3*") ribose_atoms[2] = residue_atoms[i];
            if (atm_name == " C3'") ribose_atoms[2] = residue_atoms[i];
            if (atm_name == " C4*") ribose_atoms[3] = residue_atoms[i];
            if (atm_name == " C4'") ribose_atoms[3] = residue_atoms[i];
            if (atm_name == " O4*") ribose_atoms[4] = residue_atoms[i];
            if (atm_name == " O4'") ribose_atoms[4] = residue_atoms[i];
         }
      }
      if (! (ribose_atoms[0] && ribose_atoms[1] && ribose_atoms[2] && ribose_atoms[3] && ribose_atoms[4])) {
         std::string mess = "Not all atoms found in ribose.";
         throw std::runtime_error(mess);
      } else {
         for (int i_oop_atom=0; i_oop_atom<5; i_oop_atom++) {
            clipper::Coord_orth c(ribose_atoms[i_oop_atom]->x,
                                  ribose_atoms[i_oop_atom]->y,
                                  ribose_atoms[i_oop_atom]->z);
            ribose_atoms_coords.push_back(c);
         }
         // oop: out of plane distance
         std::vector<std::pair<float, float> > pucker_distortion_and_oop_d(5);
         for (int i_oop_atom=0; i_oop_atom<5; i_oop_atom++) {
            std::vector<mmdb::Atom *> plane_atom;
            std::vector<clipper::Coord_orth> plane_atom_coords;
            for (int i=0; i<5; i++) {
               if (i != i_oop_atom) {
                  clipper::Coord_orth c(ribose_atoms[i]->x, ribose_atoms[i]->y, ribose_atoms[i]->z);
                  plane_atom.push_back(ribose_atoms[i]);
                  plane_atom_coords.push_back(c);
               }
            }
            // plane atom is now filled with 4 atoms from which the plane
            // should be calculated.
            clipper::Coord_orth pt(ribose_atoms[i_oop_atom]->x,
                                   ribose_atoms[i_oop_atom]->y,
                                   ribose_atoms[i_oop_atom]->z);
            // lsq_plane_deviation returns pair(out-of-plane-dist, rms_deviation_plane);
            std::pair<double, double> dev =
               coot::lsq_plane_deviation(plane_atom_coords, pt);
            pucker_distortion_and_oop_d[i_oop_atom] = dev;
         }

         // Find the biggest out-of-plane distance.  That is the pucker
         // of this ribose.
         puckered_atom_ = coot::pucker_analysis_info_t::NONE;
         std::pair<float, float> most_deviant(0,0);
         for (int i_oop_atom=0; i_oop_atom<5; i_oop_atom++) {
            //          std::cout << "   pucker_distortion_and_oop_d["
            //                    << i_oop_atom << "] "
            //                    << pucker_distortion_and_oop_d[i_oop_atom].first  << " "
            //                    << pucker_distortion_and_oop_d[i_oop_atom].second  << std::endl;
            if (fabs(pucker_distortion_and_oop_d[i_oop_atom].first) > fabs(most_deviant.first)) {
               most_deviant = pucker_distortion_and_oop_d[i_oop_atom];
               puckered_atom_ = possible_puckers[i_oop_atom];
            }
         }
         out_of_plane_distance = most_deviant.first;
         plane_distortion = most_deviant.second;
      }
   } else {
      // we throw on failure
      std::string mess = "base lsq plane has no value";
      throw std::runtime_error(mess);
   }
}

#include "json.hpp" // Assumes nlohmann/json.hpp is available
using json = nlohmann::json;

std::string
coot::pucker_analysis_info_t::to_json() const {

   json j;
   json j_altconf = altconf;
   json j_plane_distortion = plane_distortion;
   json j_out_of_plane_distance = out_of_plane_distance;
   json j_puckered_atom = puckered_atom();
   j["altconf"]               = j_altconf;
   j["plane_distortion"]      = j_plane_distortion;
   j["out_of_plane_distance"] = j_out_of_plane_distance;
   j["puckered_atom"]         = j_puckered_atom;

   // json["plane_distortion"] = plane_distortion;
   // json["out_of_plane_distance"] = out_of_plane_distance;
   // json["puckered_atom"] = puckered_atom();
   std::string s = j.dump(4);
   return s;
}

double coot::pucker_analysis_info_gemmi_t::compute_pucker_parameter(gemmi::Residue *residue) {
   const gemmi::Atom *o4_prime = residue->find_atom("O4'", '*');
   const gemmi::Atom *c1_prime = residue->find_atom("C1'", '*');
   const gemmi::Atom *c2_prime = residue->find_atom("C2'", '*');
   const gemmi::Atom *c3_prime = residue->find_atom("C3'", '*');
   const gemmi::Atom *c4_prime = residue->find_atom("C4'", '*');

   double v_0 = gemmi::calculate_dihedral_from_atoms(c4_prime, o4_prime, c1_prime, c2_prime);
   double v_1 = gemmi::calculate_dihedral_from_atoms(o4_prime, c1_prime, c2_prime, c3_prime);
   double v_2 = gemmi::calculate_dihedral_from_atoms(c1_prime, c2_prime, c3_prime, c4_prime);
   double v_3 = gemmi::calculate_dihedral_from_atoms(c2_prime, c3_prime, c4_prime, o4_prime);
   double v_4 = gemmi::calculate_dihedral_from_atoms(c3_prime, c4_prime, o4_prime, c1_prime);

   constexpr double deg36_r = 36.0 * M_PI / 180.0;
   constexpr double deg72_r = 72.0 * M_PI / 180.0;

   double numerator = (v_4 + v_1) - (v_3 + v_0);
   double denominator = 2.0 * v_2 * (std::sin(deg36_r) + std::sin(deg72_r));

   double P_rad = std::atan2(numerator, denominator);
   double P_deg = P_rad * 180.0 / M_PI;

   // double taum = v_2 / std::cos(P_rad);

   if (P_deg < 0)
      P_deg += 360.0;
   return P_deg;
}
coot::pucker_analysis_info_gemmi_t::pucker_analysis_info_gemmi_t(gemmi::Residue *residue, std::string altconf_in) {

   auto get_base_lsq_plane = [] (const std::vector<clipper::Coord_orth> &coords) ->
      std::optional<lsq_plane_info_t> {
      if (coords.size() < 3) return std::nullopt;
      lsq_plane_info_t lsq_plane(coords);
      return lsq_plane;
   };

   out_of_plane_distance = 0.0;
   plane_distortion = 0.0;
   N1_or_9  = nullptr;
   C1_prime = nullptr;
   puckered_atom_ = NONE;
   altconf = altconf_in;

   double P_deg = compute_pucker_parameter(residue);
   pucker_name = classify_pucker(P_deg);

   assign_base_atom_coords(residue);

   std::optional<lsq_plane_info_t> lsq_plane = get_base_lsq_plane(base_atoms_coords);
   if (!lsq_plane.has_value())
      throw std::runtime_error("base lsq plane has no value");

   markup_info.base_ring_centre = lsq_plane.value().centre();
   markup_info.base_ring_normal = lsq_plane.value().normal();

   char altloc_char = altconf_in.empty() ? '\0' : altconf_in[0];

   // find the phosphorus atom in this residue
   for (const gemmi::Atom &atom : residue->atoms) {
      if (atom.name == "P" && atom.altloc == altloc_char) {
         clipper::Coord_orth p(atom.pos.x, atom.pos.y, atom.pos.z);
         markup_info.phosphorus_position = p;
         markup_info.projected_point = lsq_plane.value().projected_point(p);
      }
   }

   // find the ribose atoms
   const std::array<std::string, 5> ribose_names = {"C1'", "C2'", "C3'", "C4'", "O4'"};
   std::vector<const gemmi::Atom *> ribose_atoms(5, nullptr);
   for (const gemmi::Atom &atom : residue->atoms) {
      if (atom.altloc != altloc_char) continue;
      for (int i = 0; i < 5; i++) {
         if (atom.name == ribose_names[i])
            ribose_atoms[i] = &atom;
      }
   }

   for (int i = 0; i < 5; i++) {
      if (!ribose_atoms[i])
         throw std::runtime_error("Not all atoms found in ribose.");
      const gemmi::Atom &a = *ribose_atoms[i];
      ribose_atoms_coords.emplace_back(a.pos.x, a.pos.y, a.pos.z);
   }

   const std::array<PUCKERED_ATOM_T, 5> possible_puckers = {
      C1_PRIME, C2_PRIME, C3_PRIME, C4_PRIME, O4_PRIME
   };
   std::pair<float, float> most_deviant(0, 0);
   for (int i_oop = 0; i_oop < 5; i_oop++) {
      std::vector<clipper::Coord_orth> plane_coords;
      for (int i = 0; i < 5; i++) {
         if (i != i_oop)
            plane_coords.push_back(ribose_atoms_coords[i]);
      }
      clipper::Coord_orth pt = ribose_atoms_coords[i_oop];
      std::pair<double, double> dev = coot::lsq_plane_deviation(plane_coords, pt);
      if (std::fabs(dev.first) > std::fabs(most_deviant.first)) {
         most_deviant = {static_cast<float>(dev.first), static_cast<float>(dev.second)};
         puckered_atom_ = possible_puckers[i_oop];
      }
   }
   out_of_plane_distance = most_deviant.first;
   plane_distortion = most_deviant.second;
}

float
coot::pucker_analysis_info_gemmi_t::phosphate_distance_to_base_plane(gemmi::Residue *following_res) {

   float oop = 0.0;
   bool found = false;
   char altloc_char = altconf.empty() ? '\0' : altconf[0];

   for (const gemmi::Atom &atom : following_res->atoms) {
      if (atom.name == "P" && atom.altloc == altloc_char) {
         if (base_atoms_coords.size() < 4) {
            std::string m = "Failed to find base atoms. Found ";
            m += std::to_string(base_atoms_coords.size());
            m += " atoms.";
            throw std::runtime_error(m);
         }
         clipper::Coord_orth pt(atom.pos.x, atom.pos.y, atom.pos.z);
         std::pair<double, double> oop_plus_dev = coot::lsq_plane_deviation(base_atoms_coords, pt);
         oop = static_cast<float>(oop_plus_dev.first);
         markup_info.phosphorus_position = pt;
         lsq_plane_info_t lsq_plane(base_atoms_coords);
         markup_info.projected_point = lsq_plane.projected_point(pt);
         found = true;
         break;
      }
   }
   if (!found)
      throw std::runtime_error("Failed to find following phosphate");
   return oop;
}

float
coot::pucker_analysis_info_gemmi_t::phosphate_distance(gemmi::Residue *following_res) {

   if (!C1_prime)
      throw std::runtime_error("C1' not found in this residue");
   if (!N1_or_9)
      throw std::runtime_error("N1/N9 not found in this residue");

   float d = 0.0;
   bool found = false;
   char altloc_char = altconf.empty() ? '\0' : altconf[0];

   for (const gemmi::Atom &atom : following_res->atoms) {
      if (atom.name == "P" && atom.altloc == altloc_char) {
         clipper::Coord_orth P_pt(atom.pos.x, atom.pos.y, atom.pos.z);
         clipper::Coord_orth N_pt(N1_or_9->pos.x,  N1_or_9->pos.y,  N1_or_9->pos.z);
         clipper::Coord_orth C_pt(C1_prime->pos.x, C1_prime->pos.y, C1_prime->pos.z);
         clipper::Coord_orth CN = N_pt - C_pt;
         clipper::Coord_orth CP = P_pt - C_pt;
         double CN_d = clipper::Coord_orth::length(N_pt, C_pt);
         double CP_d = clipper::Coord_orth::length(P_pt, C_pt);
         if (CN_d > 0.0 && CP_d > 0.0) {
            double cos_alpha = clipper::Coord_orth::dot(CN, CP) / (CN_d * CP_d);
            double alpha = std::acos(cos_alpha);
            d = static_cast<float>(std::sin(M_PI - alpha) * CP_d);
            found = true;
         }
         break;
      }
   }
   if (!found)
      throw std::runtime_error("P not found in following residue");
   return d;
}

std::string
coot::pucker_analysis_info_gemmi_t::puckered_atom() const {

   switch (puckered_atom_) {
      case C1_PRIME: return "C1'";
      case C2_PRIME: return "C2'";
      case C3_PRIME: return "C3'";
      case C4_PRIME: return "C4'";
      case O4_PRIME: return "O4'";
      default:       return "----";
   }
}

json coot::pucker_analysis_info_gemmi_t::to_json() const {

   json j;
   j["altconf"]               = altconf;
   j["pucker_name"]           = pucker_name;
   j["plane_distortion"]      = plane_distortion;
   j["out_of_plane_distance"] = out_of_plane_distance;
   j["puckered_atom"]         = puckered_atom();
   return j;
}
std::string coot::pucker_analysis_info_gemmi_t::classify_pucker(double P_deg) {
   const std::array<std::string, 10> puckers = {
         "C3'-endo",
         "C4'-exo",
         "O4'-endo",
         "C1'-exo",
         "C2'-endo",
         "C3'-exo",
         "C4'-endo",
         "O4'-exo",
         "C1'-endo",
         "C2'-exo"
     };

   P_deg = std::fmod(P_deg, 360.0);
   if (P_deg < 0.0)
      P_deg += 360.0;

   const auto index =static_cast<std::size_t>(std::floor((P_deg + 18.0) / 36.0)) % puckers.size();
   return puckers[index];
}

void
coot::pucker_analysis_info_gemmi_t::assign_base_atom_coords(gemmi::Residue *residue) {

   // gemmi atom names are unpadded (e.g. "N1", "C2", "O4'") unlike the 4-char PDB-padded mmdb names
   const std::vector<std::string> cytidine_base_names  = {"N1","C2","N3","C4","C5","C6","O2","N4"};
   const std::vector<std::string> uracil_base_names    = {"N1","C2","N3","C4","C5","C6","O2","O4"};
   const std::vector<std::string> adenine_base_names   = {"N9","C8","N7","C5","C4","N1","C2","N3","C6","N6"};
   const std::vector<std::string> guanine_base_names   = {"N9","C8","N7","C5","C4","N1","C2","N3","C6","O6","N2"};
   const std::vector<std::string> thymine_base_names   = {"N1","C2","N3","C4","C5","C6","O2","O4","C7"};

   std::vector<std::string> base_names;
   const std::string &res_name = residue->name;
   if (res_name == "C" || res_name == "Cr" || res_name == "CYT" || res_name == "DC")
      base_names = cytidine_base_names;
   else if (res_name == "U" || res_name == "Ur" || res_name == "URA")
      base_names = uracil_base_names;
   else if (res_name == "A" || res_name == "Ar" || res_name == "ADE" || res_name == "DA")
      base_names = adenine_base_names;
   else if (res_name == "G" || res_name == "Gr" || res_name == "GUA" || res_name == "DG")
      base_names = guanine_base_names;
   else if (res_name == "T" || res_name == "Tr" || res_name == "THY" || res_name == "DT")
      base_names = thymine_base_names;

   char altloc_char = altconf.empty() ? '\0' : altconf[0];

   for (const gemmi::Atom &atom : residue->atoms) {
      if (atom.altloc != altloc_char) continue;
      if (atom.name == "N1" || atom.name == "N9")
         N1_or_9 = &atom;
      if (atom.name == "C1'")
         C1_prime = &atom;
      if (!base_names.empty()) {
         for (const std::string &bn : base_names) {
            if (atom.name == bn) {
               base_atoms_coords.emplace_back(atom.pos.x, atom.pos.y, atom.pos.z);
               break;
            }
         }
      }
   }
}

void
coot::pucker_analysis_info_t::assign_base_atom_coords(mmdb::Residue *residue_p) {

   std::vector<std::string> cytidine_base_names;
   std::vector<std::string> uracil_base_names;
   std::vector<std::string> adenine_base_names;
   std::vector<std::string> guanine_base_names;
   std::vector<std::string> thymine_base_names;

   cytidine_base_names.push_back(" N1 ");
   cytidine_base_names.push_back(" C2 ");
   cytidine_base_names.push_back(" N3 ");
   cytidine_base_names.push_back(" C4 ");
   cytidine_base_names.push_back(" C5 ");
   cytidine_base_names.push_back(" C6 ");
   cytidine_base_names.push_back(" O2 ");
   cytidine_base_names.push_back(" N4 ");

   uracil_base_names.push_back(" N1 ");
   uracil_base_names.push_back(" C2 ");
   uracil_base_names.push_back(" N3 ");
   uracil_base_names.push_back(" C4 ");
   uracil_base_names.push_back(" C5 ");
   uracil_base_names.push_back(" C6 ");
   uracil_base_names.push_back(" O2 ");
   uracil_base_names.push_back(" O4 ");

   adenine_base_names.push_back(" N9 ");
   adenine_base_names.push_back(" C8 ");
   adenine_base_names.push_back(" N7 ");
   adenine_base_names.push_back(" C5 ");
   adenine_base_names.push_back(" C4 ");
   adenine_base_names.push_back(" N1 ");
   adenine_base_names.push_back(" C2 ");
   adenine_base_names.push_back(" N3 ");
   adenine_base_names.push_back(" C6 ");
   adenine_base_names.push_back(" N6 ");

   guanine_base_names.push_back(" N9 ");
   guanine_base_names.push_back(" C8 ");
   guanine_base_names.push_back(" N7 ");
   guanine_base_names.push_back(" C5 ");
   guanine_base_names.push_back(" C4 ");
   guanine_base_names.push_back(" N1 ");
   guanine_base_names.push_back(" C2 ");
   guanine_base_names.push_back(" N3 ");
   guanine_base_names.push_back(" C6 ");
   guanine_base_names.push_back(" O6 ");
   guanine_base_names.push_back(" N2 ");

   thymine_base_names.push_back(" N1 ");
   thymine_base_names.push_back(" C2 ");
   thymine_base_names.push_back(" N3 ");
   thymine_base_names.push_back(" C4 ");
   thymine_base_names.push_back(" C5 ");
   thymine_base_names.push_back(" C6 ");
   thymine_base_names.push_back(" O2 ");
   thymine_base_names.push_back(" O4 ");
   thymine_base_names.push_back(" C5M");


   mmdb::PPAtom residue_atoms = NULL;
   int n_residue_atoms;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);

   // Assign N1_or_9 and C1_prime
   for (int i=0; i<n_residue_atoms; i++) {
      std::string atom_name(residue_atoms[i]->name);
      std::string alt_name(residue_atoms[i]->altLoc);
      if (alt_name == altconf) {
         if (atom_name == " N1 ")
            N1_or_9 = residue_atoms[i];
         if (atom_name == " N9 ")
            N1_or_9 = residue_atoms[i];
         if (atom_name == " C1*")
            C1_prime = residue_atoms[i];
         if (atom_name == " C1'")
            C1_prime = residue_atoms[i];
      }
   }

   // Fill base_names according to residue type/name.  If base_name is
   // empty after setting, just fall out (an exception is thrown in
   // the constructor if there are not enough base name atoms.

   std::vector<std::string> base_names;

   std::string residue_name(residue_p->GetResName());

   // current names
   if (residue_name == "C") base_names = cytidine_base_names;
   if (residue_name == "U") base_names = uracil_base_names;
   if (residue_name == "A") base_names = adenine_base_names;
   if (residue_name == "G") base_names = guanine_base_names;
   // old names
   if (residue_name == "Cr") base_names = cytidine_base_names;
   if (residue_name == "Ur") base_names = uracil_base_names;
   if (residue_name == "Ar") base_names = adenine_base_names;
   if (residue_name == "Gr") base_names = guanine_base_names;
   // modern (3.x) RNA base names
   if (residue_name == "CYT") base_names = cytidine_base_names;
   if (residue_name == "URA") base_names = uracil_base_names;
   if (residue_name == "ADE") base_names = adenine_base_names;
   if (residue_name == "GUA") base_names = guanine_base_names;

   if (base_names.size() > 0) {
      for (int i=0; i<n_residue_atoms; i++) {
         std::string atm_name(residue_atoms[i]->name);
         std::string alt_name(residue_atoms[i]->altLoc);
         for (unsigned int j=0; j<base_names.size(); j++) {
            if (base_names[j] == atm_name) {
               base_atoms_coords.push_back(clipper::Coord_orth(residue_atoms[i]->x,
                                                               residue_atoms[i]->y,
                                                               residue_atoms[i]->z));
            }
         }
      }
   }
} 


// Use the 3' phosphate of the following residue to calculate its out
// of plane distance (the plane being the base plane).  Decide from
// that if this should have been 3' or 2'.  Check vs the actual
// puckering.
//
// Throw an exception if we can't do this.
// 
float 
coot::pucker_analysis_info_t::phosphate_distance_to_base_plane(mmdb::Residue *following_res) {

   float oop = 0.0;
   mmdb::PPAtom residue_atoms = NULL;
   int n_residue_atoms;
   bool found = 0;

   following_res->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int i=0; i<n_residue_atoms; i++) {
      std::string atm_name(residue_atoms[i]->name);
      std::string alt_name(residue_atoms[i]->altLoc);
      if (atm_name == " P  ") { 
         if (altconf == alt_name) {
            clipper::Coord_orth pt(residue_atoms[i]->x,
                                   residue_atoms[i]->y,
                                   residue_atoms[i]->z);
            // lsq_plane_deviation returns pair(out-of-plane-dist, rms_deviation_plane);

            if (base_atoms_coords.size() < 4) {

               // construct an error message and throw an exception.
               // 
               std::string m = "Failed to find base atoms. Found ";
               m += coot::util::int_to_string(base_atoms_coords.size());
               m += " atoms. ";
               throw std::runtime_error(m);
            } else {
               std::pair<double, double> oop_plus_dev =
                  coot::lsq_plane_deviation(base_atoms_coords, pt);
               oop = oop_plus_dev.first;
               // overwrite the markup info so that phosphorus_position and
               // projected_point refer to the following residue's P (which is
               // the P used for the pukka-pucker test), not the current
               // residue's own P.
               markup_info.phosphorus_position = pt;
               lsq_plane_info_t lsq_plane(base_atoms_coords);
               markup_info.projected_point = lsq_plane.projected_point(pt);
               found = 1;
               break;
            }
         }
      }
   }
   if (found == 0) {
      throw std::runtime_error("Failed to find following phosphate");
   }
   return oop;
}

// Throw an exception if the reference atoms are not found.
float 
coot::pucker_analysis_info_t::phosphate_distance(mmdb::Residue *following_res) {

   if (! C1_prime) { 
      std::string mess = "C1*/C1' not found in this residue";
      throw std::runtime_error(mess);
   }
   if (! N1_or_9) { 
      std::string mess = "N1/N9 not found in this residue";
      throw std::runtime_error(mess);
   }

   //                     X             
   //                    / \  90 degrees       ;
   //                 d / X \                  ;
   //                  /     \                 ;
   //                 /       \                       ;
   //                /         \               ; 
   //               /   pi-alpha\              ;
   //              --------------\C1'               ;
   //            P           alpha\            ; 
   //                              \           ;
   //                               \          ; 
   //                                \         ;
   //                                 \               ;
   //                                  \       ; 
   //                                   \ N1   ;
   
   
   float d = 0.0;
   mmdb::PPAtom residue_atoms = NULL;
   int n_residue_atoms;
   bool found = 0;
   following_res->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int i=0; i<n_residue_atoms; i++) {
      std::string atm_name(residue_atoms[i]->name);
      std::string alt_name(residue_atoms[i]->altLoc);
      if (atm_name == " P  ") { 
         if (altconf == alt_name) {
            clipper::Coord_orth P_pt(residue_atoms[i]->x,
                                     residue_atoms[i]->y,
                                     residue_atoms[i]->z);
            clipper::Coord_orth N_pt( N1_or_9->x,  N1_or_9->y,  N1_or_9->z);
            clipper::Coord_orth C_pt(C1_prime->x, C1_prime->y, C1_prime->z);
            clipper::Coord_orth CN = N_pt - C_pt;
            clipper::Coord_orth CP = P_pt - C_pt;

            double CN_d = clipper::Coord_orth::length(N_pt, C_pt);
            double CP_d = clipper::Coord_orth::length(P_pt, C_pt);

            if (CN_d > 0.0) { 
               if (CP_d > 0.0) { 
                  found = 1;
                  double cos_alpha = clipper::Coord_orth::dot(CN, CP)/(CN_d*CP_d);
                  
                  double alpha = acos(cos_alpha);
                  double sin_pi_minus_alpha = sin(M_PI - alpha);
                  d = sin_pi_minus_alpha * CP_d;
               }
            }
         }
      }
   }

   if (! found) {
      std::string mess = "P not found in this residue";
      throw std::runtime_error(mess);
   }
   return d;
} 


std::string
coot::pucker_analysis_info_t::puckered_atom() const {

   std::string s;
   if (puckered_atom_ == coot::pucker_analysis_info_t::C2_PRIME)
      s = " C2'";
   if (puckered_atom_ == coot::pucker_analysis_info_t::C3_PRIME)
      s = " C3'";
   if (puckered_atom_ == coot::pucker_analysis_info_t::NONE)
      s = "----";
   if (puckered_atom_ == coot::pucker_analysis_info_t::C1_PRIME)
      s = " C1'";
   if (puckered_atom_ == coot::pucker_analysis_info_t::C4_PRIME)
      s = " C4'";
   if (puckered_atom_ == coot::pucker_analysis_info_t::O4_PRIME)
      s = " O4'";
   return s;
} 


// return "" on no canonical name found
std::string
coot::util::canonical_base_name(const std::string &res_name_in, base_t rna_or_dna) {

   if (rna_or_dna == coot::RNA) {
      if (res_name_in == "C")
         return "C";
      if (res_name_in == "A")
         return "A";
      if (res_name_in == "G")
         return "G";
      if (res_name_in == "T")
         return "T";
      if (res_name_in == "U")
         return "U";
      if (res_name_in == "Cr")
         return "C";
      if (res_name_in == "Ar")
         return "A";
      if (res_name_in == "Gr")
         return "G";
      if (res_name_in == "Tr")
         return "T";
      if (res_name_in == "Ur")
         return "U";
      if (res_name_in == "Cd")
         return "C";
      if (res_name_in == "Ad")
         return "A";
      if (res_name_in == "Gd")
         return "G";
      if (res_name_in == "Td")
         return "T";
      if (res_name_in == "Ud")
         return "U";
   }
   
   if (rna_or_dna == coot::DNA) {
      if (res_name_in == "C")
         return "DC";
      if (res_name_in == "A")
         return "DA";
      if (res_name_in == "G")
         return "DG";
      if (res_name_in == "T")
         return "DT";
      if (res_name_in == "U")
         return "DU";
      if (res_name_in == "Cd")
         return "DC";
      if (res_name_in == "Ad")
         return "DA";
      if (res_name_in == "Gd")
         return "DG";
      if (res_name_in == "Td")
         return "DT";
      if (res_name_in == "Ud")
         return "DU";
      if (res_name_in == "Cr")
         return "DC";
      if (res_name_in == "Ar")
         return "DA";
      if (res_name_in == "Gr")
         return "DG";
      if (res_name_in == "Tr")
         return "DT";
      if (res_name_in == "Ur")
         return "DU";
   }

   return "";
}
