#pragma once

#include <gemmi/model.hpp>
#include <utility>
#include <string>
#include "clipper/core/coords.h"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__is_nucleotide_by_dict/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"
#include "geometry/protein-geometry.hh"

// Standalone copy of coot::restraints_container_t::link_restraints_counts
// (which is private in the original class) so the _gemmi function can return it.
struct link_restraints_counts_gemmi {
   link_restraints_counts_gemmi() {
      init();
   }
   void init() {
      n_link_bond_restr = 0;
      n_link_angle_restr = 0;
      n_link_trans_peptide = 0;
      n_link_torsion_restr = 0;
      n_link_plane_restr = 0;
      n_link_improper_dihedral_restr = 0;
      link_type = "link";
   }
   std::string link_type;
   unsigned int n_link_bond_restr;
   unsigned int n_link_angle_restr;
   unsigned int n_link_plane_restr;
   unsigned int n_link_torsion_restr;
   unsigned int n_link_trans_peptide;
   unsigned int n_link_improper_dihedral_restr;
};

inline std::pair<bool, link_restraints_counts_gemmi>
try_make_phosphodiester_link_ng_gemmi(
    const coot::protein_geometry &geom,
    std::pair<bool, const gemmi::Residue *> res_1_pair,
    std::pair<bool, const gemmi::Residue *> res_2_pair) {

   bool use_distance_cut_off = false;
   const float distance_cut_off = 5.0;
   const float distance_cut_off_srd = distance_cut_off * distance_cut_off;

   const gemmi::Residue *res_1 = res_1_pair.second;
   const gemmi::Residue *res_2 = res_2_pair.second;

   bool status = false;
   link_restraints_counts_gemmi lrc;

   if (coot::util::is_nucleotide_by_dict_gemmi(*res_1, res_1->name)) {
      if (coot::util::is_nucleotide_by_dict_gemmi(*res_2, res_2->name)) {

         for (const gemmi::Atom &at_1 : res_1->atoms) {
            std::string at_name_1 = at_1.name;
            if (at_name_1 == " O3'") {
               char alt_conf_1 = at_1.altloc;
               for (const gemmi::Atom &at_2 : res_2->atoms) {
                  std::string at_name_2 = at_2.name;
                  if (at_name_2 == " P  ") {
                     char alt_conf_2 = at_2.altloc;
                     if (alt_conf_1 == alt_conf_2 || alt_conf_1 == '\0' || alt_conf_2 == '\0') {

                        if (use_distance_cut_off) {
                           int res_no_1 = res_1->seqid.num.value;
                           int res_no_2 = res_2->seqid.num.value;
                           if ((res_no_2 - res_no_1) > 1) {
                              clipper::Coord_orth pt_1 = coot::co_gemmi(&at_1);
                              clipper::Coord_orth pt_2 = coot::co_gemmi(&at_2);
                              if ((pt_2 - pt_1).lengthsq() > distance_cut_off_srd)
                                 continue;
                           }
                        }

                        std::string link_type = "p";
                        bool is_fixed_first_residue  = res_1_pair.first;
                        bool is_fixed_second_residue = res_2_pair.second;

                        // make_link_restraints_for_link_ng has no gemmi port;
                        // original just adds restraints and returns counts.
                        // We set status=true; lrc stays zeroed.
                        status = true;
                     }
                  }
               }
            }
         }
      }
   }

   return std::pair<bool, link_restraints_counts_gemmi>(status, lrc);
}