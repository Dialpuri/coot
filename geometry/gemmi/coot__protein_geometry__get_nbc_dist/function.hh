#pragma once

#include <gemmi/elem.hpp>
#include <coot/geometry/protein-geometry.hh>
#include <coot/geometry/hb-types.hh>
#include <coot/geometry/energy-lib.hh>
#include <utility>
#include <string>

namespace coot {

inline std::pair<bool, double>
get_nbc_dist_gemmi(const coot::protein_geometry& geom,
                   const std::string& energy_type_1,
                   const std::string& energy_type_2,
                   bool in_same_residue_flag,
                   bool in_same_ring_flag) {

   std::pair<bool, double> r(false, 0.0);

   coot::energy_lib_atom atom_1 = geom.get_energy_lib_atom(energy_type_1);
   coot::energy_lib_atom atom_2 = geom.get_energy_lib_atom(energy_type_2);

   if (atom_1.type.empty() || atom_2.type.empty()) {
      return r;
   }

   r.first = true;

   float radius_1;
   float radius_2;

   gemmi::Element elem_1(atom_1.element);
   if (elem_1.is_metal()) {
      radius_1 = atom_1.ion_radius;
   } else {
      radius_1 = atom_1.vdw_radius;
   }

   gemmi::Element elem_2(atom_2.element);
   if (elem_2.is_metal()) {
      radius_2 = atom_2.ion_radius;
   } else {
      radius_2 = atom_2.vdw_radius;
   }

   r.second = radius_1 + radius_2;

   if (in_same_residue_flag) {
      r.second *= 0.84;
   }

   if (in_same_ring_flag) {
      auto is_ring_atom = [](const std::string& name) {
         return name == "CR15" || name == "CR16" || name == "CR1"  ||
                name == "CR6"  || name == "CR5"  || name == "CR56" ||
                name == "CR66" || name == "NPA"  || name == "NPB"  ||
                name == "NRD5" || name == "NRD6" || name == "NR15" ||
                name == "NR16" || name == "NR6"  || name == "NR5";
      };
      if (is_ring_atom(energy_type_1) && is_ring_atom(energy_type_2)) {
         r.second = 2.2;
      }
   }

   // Hydrogen bonds can be closer
   if ((atom_1.hb_type == HB_DONOR ||
        atom_1.hb_type == HB_BOTH  ||
        atom_1.hb_type == HB_HYDROGEN) &&
       (atom_2.hb_type == HB_ACCEPTOR ||
        atom_2.hb_type == HB_BOTH)) {
      r.second -= 0.5;
      if (atom_1.hb_type == HB_HYDROGEN)
         r.second -= 0.3;
   }

   if ((atom_2.hb_type == HB_DONOR ||
        atom_2.hb_type == HB_BOTH  ||
        atom_2.hb_type == HB_HYDROGEN) &&
       (atom_1.hb_type == HB_ACCEPTOR ||
        atom_1.hb_type == HB_BOTH)) {
      r.second -= 0.5;
      if (atom_2.hb_type == HB_HYDROGEN)
         r.second -= 0.3;
   }

   return r;
}

} // namespace coot