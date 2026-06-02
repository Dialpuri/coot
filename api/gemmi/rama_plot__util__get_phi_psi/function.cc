#include "function.hh"
#include "clipper/clipper.h"
#include <cstdio>

std::pair<bool, rama_plot::phi_psi_t>
rama_plot::util::get_phi_psi_gemmi(gemmi::CRA residue_0, gemmi::CRA residue_1, gemmi::CRA residue_2) {

   bool is_valid_flag = false;
   bool is_pre_pro = false;
   rama_plot::phi_psi_t phi_psi; // part of the returned value
   int natom = 0;
   int ires = 0;
   clipper::Coord_orth c_prev, n_this, ca_this, c_this, n_next;

   // Extract C atom from residue_0
   if (residue_0.residue) {
      for (const auto& atom : residue_0.residue->atoms) {
         std::string atom_name = atom.name;
         if (atom_name == "C") {
            c_prev = clipper::Coord_orth(atom.pos.x, atom.pos.y, atom.pos.z);
            natom++;
         }
      }
   }

   std::string res_name_1;
   if (residue_1.residue) {
      res_name_1 = residue_1.residue->name;
      ires = residue_1.residue->seqid.num.value;
      for (const auto& atom : residue_1.residue->atoms) {
         std::string atom_name = atom.name;
         if (atom_name == "C") {
            c_this = clipper::Coord_orth(atom.pos.x, atom.pos.y, atom.pos.z);
            natom++;
         }
         if (atom_name == "CA") {
            ca_this = clipper::Coord_orth(atom.pos.x, atom.pos.y, atom.pos.z);
            natom++;
         }
         if (atom_name == "N") {
            n_this = clipper::Coord_orth(atom.pos.x, atom.pos.y, atom.pos.z);
            natom++;
         }
      }
   }

   // Check if residue_2 is PRO
   if (residue_2.residue && std::string(residue_2.residue->name) == "PRO")
      is_pre_pro = true;

   if (residue_2.residue) {
      for (const auto& atom : residue_2.residue->atoms) {
         std::string atom_name = atom.name;
         if (atom_name == "N") {
            n_next = clipper::Coord_orth(atom.pos.x, atom.pos.y, atom.pos.z);
            natom++;
         }
      }
   }

   if (natom == 5 && residue_1.residue && residue_1.chain) {
      char num[30];
      std::snprintf(num, 20, "%d", ires);
      std::string label(num);
      std::string segid = residue_1.chain->name;
      
      // Get insertion code - normalize: MMDB uses "" for no icode, gemmi uses ' '
      std::string inscode;
      if (residue_1.residue->seqid.icode == ' ') {
         inscode = "";
      } else {
         inscode = std::string(1, residue_1.residue->seqid.icode);
      }
      
      label += inscode;
      label += " ";
      label += segid;
      label += " ";
      label += residue_1.residue->name;

      double phi = clipper::Util::rad2d(ca_this.torsion(c_prev, n_this, ca_this, c_this));
      double psi = clipper::Util::rad2d(ca_this.torsion(n_this, ca_this, c_this, n_next));

      phi_psi = rama_plot::phi_psi_t(phi, psi,
                                     residue_1.residue->name,
                                     label.c_str(),
                                     ires,
                                     inscode,
                                     segid,
                                     is_pre_pro);
      // peptide bonding atoms have to be within 2.0A, or this is not
      // a valid peptide.
      //
      double dist_1 = clipper::Coord_orth::length(c_prev, n_this);
      double dist_2 = clipper::Coord_orth::length(c_this, n_next);

      if (dist_1 < 2.0)
         if (dist_2 < 2.0)
            is_valid_flag = true;

   }

   rama_plot::phi_psi_t phi_psi_with_residues(phi_psi);
   phi_psi_with_residues.type = clipper::Ramachandran::All;
   if (res_name_1 == "GLY") {
      phi_psi_with_residues.type = clipper::Ramachandran::Gly2;
   } else {
      if (res_name_1 == "PRO") {
         phi_psi_with_residues.type = clipper::Ramachandran::Pro2;
      } else {
         if (res_name_1 == "ILE" || res_name_1 == "VAL") {
            phi_psi_with_residues.type = clipper::Ramachandran::IleVal2;
         } else {
            phi_psi_with_residues.type = clipper::Ramachandran::NoGPIVpreP2;
         }
      }
   }
   if (is_pre_pro) {
      phi_psi_with_residues.type = clipper::Ramachandran::PrePro2;
   }

   return std::pair<bool, rama_plot::phi_psi_t>(is_valid_flag, phi_psi_with_residues);

}