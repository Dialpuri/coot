#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

#include "geometry/protein-geometry.hh"

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__is_linked/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__place_hydrogen_by_connected_atom_energy_type__593ada/gemmi/function.hh"

namespace coot::reduce {

void hydrogen_placement_by_dictionary_gemmi(
    const dictionary_residue_restraints_t &rest,
    gemmi::CRA cra,
    const gemmi::Structure &structure,
    double bl_aliph,
    double bl_arom,
    double bl_amino,
    double bl_oh,
    double bl_sh)
{
   std::vector<std::string> done_atom_name_list;

   for (unsigned int iat = 0; iat < rest.atom_info.size(); iat++) {
      if (rest.atom_info[iat].is_hydrogen()) {
         const std::string &H_at_name = rest.atom_info[iat].atom_id_4c;
         // if we haven't done it already...
         if (std::find(done_atom_name_list.begin(), done_atom_name_list.end(), H_at_name) == done_atom_name_list.end()) {
            // skip the HO3' on RNA and DNA
            if ((rest.residue_info.group == "DNA" || rest.residue_info.group == "RNA") &&
                H_at_name == "HO3'") {
               continue;
            } else {
               // to which atom is this hydrogen connected?
               std::vector<unsigned int> neighbs = rest.neighbours(iat, false);
               if (neighbs.size() == 1) {
                  const unsigned int &iat_neighb = neighbs[0];
                  const std::string &energy_type = rest.atom_info[iat_neighb].type_energy;
                  const std::string &first_neigh = rest.atom_info[iat_neighb].atom_id_4c;
                  if (!is_linked_gemmi(first_neigh, cra, structure)) {
                     if (!energy_type.empty()) {
                        std::vector<std::string> v =
                           place_hydrogen_by_connected_atom_energy_type_gemmi(
                               energy_type, iat, iat_neighb, rest, cra.residue,
                               bl_aliph, bl_arom, bl_amino, bl_oh, bl_sh);
                        done_atom_name_list.insert(done_atom_name_list.end(), v.begin(), v.end());
                     } else {
                        // No energy type — use 2nd neighbours
                        // (original: place_hydrogen_by_connected_2nd_neighbours)
                        const std::string &H_at_name_2nd = rest.atom_info[iat].atom_id_4c;
                        const std::string &at_name_1_2nd = rest.atom_info[iat_neighb].atom_id_4c;
                        std::vector<unsigned int> iat_neighb_2 = rest.neighbours(iat_neighb, false);
                        std::string at_name_2_2nd = "";
                        if (!iat_neighb_2.empty()) {
                           at_name_2_2nd = rest.atom_info[iat_neighb_2[0]].atom_id_4c;
                        }
                        // Use add_OH_H_gemmi for 2nd-neighbour fallback
                        // (angle PI/2, torsion 0 as defaults when no energy type)
                        coot::reduce::add_OH_H_gemmi(
                            H_at_name_2nd, at_name_1_2nd, at_name_2_2nd, "",
                            bl_aliph, M_PI_2, 0.0, *cra.residue);
                        done_atom_name_list.push_back(H_at_name_2nd);
                     }
                  }
               }
            }
         }
      }
   }
}

} // namespace coot::reduce