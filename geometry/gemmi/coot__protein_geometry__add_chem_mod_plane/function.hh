#pragma once

#include <gemmi/cif.hpp>
#include <string>
#include <iostream>
#include <vector>
#include "coot/geometry/protein-geometry.hh"

namespace coot {

inline void add_chem_mod_plane_gemmi(protein_geometry& geom, const gemmi::cif::Loop& mmCIFLoop) {

   if (mmCIFLoop.tags.empty()) return;

   unsigned num_rows = mmCIFLoop.length();
   if (num_rows == 0) return;

   // Find column indices for each field (short names as in original MMDB code)
   // The original MMDB code uses short names like "mod_id" without category prefix.
   // mmdb::mmcif::Loop::GetString returns ierr=-25 when the short name doesn't match
   // a column. We replicate this: find_tag with short names returns -1 when the column
   // doesn't exist.
   int col_mod_id = mmCIFLoop.find_tag("mod_id");
   int col_plane_id = mmCIFLoop.find_tag("plane_id");
   int col_function = mmCIFLoop.find_tag("function");
   int col_atom_id = mmCIFLoop.find_tag("atom_id");
   int col_new_dist_esd = mmCIFLoop.find_tag("new_dist_esd");

   for (unsigned j = 0; j < num_rows; j++) {
      int ierr_tot = 0;
      
      std::string mod_id;
      std::string plane_id;
      std::string function;
      std::string atom_id;
      double new_dist_esd = 0.0;

      // Extract mod_id
      if (col_mod_id >= 0) {
         mod_id = mmCIFLoop.val(j, col_mod_id);
         if (mod_id == "?") {
            mod_id = "";
            ierr_tot += 1;
         }
      } else {
         ierr_tot += 1;
      }

      // Extract plane_id
      if (col_plane_id >= 0) {
         plane_id = mmCIFLoop.val(j, col_plane_id);
         if (plane_id == "?") {
            plane_id = "";
            ierr_tot += 1;
         }
      } else {
         ierr_tot += 1;
      }

      // Extract function
      if (col_function >= 0) {
         function = mmCIFLoop.val(j, col_function);
         if (function == "?") {
            function = "";
            ierr_tot += 1;
         }
      } else {
         ierr_tot += 1;
      }

      // Extract atom_id
      if (col_atom_id >= 0) {
         atom_id = mmCIFLoop.val(j, col_atom_id);
         if (atom_id == "?") {
            atom_id = "";
            ierr_tot += 1;
         }
      } else {
         ierr_tot += 1;
      }

      // Extract new_dist_esd
      if (col_new_dist_esd >= 0) {
         std::string val = mmCIFLoop.val(j, col_new_dist_esd);
         if (val != "?" && !val.empty()) {
            try {
               new_dist_esd = std::stod(val);
            } catch (...) {
               new_dist_esd = 0.0;
               ierr_tot += 1;
            }
         } else {
            new_dist_esd = 0.0;
            ierr_tot += 1;
         }
      } else {
         new_dist_esd = 0.0;
         ierr_tot += 1;
      }

      if (ierr_tot == 0 || function == "delete") {
         std::string atom_name = coot::atom_id_mmdb_expand(atom_id);
         geom.mods[mod_id].add_plane_atom(plane_id, function, atom_name, new_dist_esd);
      } else {
         // original: std::cout << "oops in add_chem_mod_plane ierr_tot is " << ierr_tot << std::endl;
         (void)ierr_tot; // suppress unused warning
      }
   }
}

} // namespace coot