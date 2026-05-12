#include "function.hh"

namespace coot {

atom_overlaps_container_t::h_bond_info_t::h_bond_info_t(int hb_1, int hb_2,
                                                        const std::string& resname_1,
                                                        const std::string& resname_2) {
   is_h_bond_H_and_acceptor = false;
   is_h_bond_donor_and_acceptor = false;
   H_is_first_atom_flag = false;
   H_is_second_atom_flag = false;
   donor_is_second_atom_flag = false;

   if (hb_1 == HB_HYDROGEN) {
      if (hb_2 == HB_ACCEPTOR || hb_2 == HB_BOTH) {
         is_h_bond_H_and_acceptor = true;
         H_is_first_atom_flag = true;
      }
   }

   if (hb_1 == HB_ACCEPTOR || hb_1 == HB_BOTH) {
      if (hb_2 == HB_HYDROGEN) {
         is_h_bond_H_and_acceptor = true;
         H_is_second_atom_flag = true;
      }
   }

   if (hb_1 == HB_DONOR || hb_1 == HB_BOTH) {
      if (hb_2 == HB_ACCEPTOR || hb_2 == HB_BOTH) {
         is_h_bond_donor_and_acceptor = true;
      }
   }

   if (hb_1 == HB_ACCEPTOR || hb_1 == HB_BOTH) {
      if (hb_2 == HB_DONOR || hb_2 == HB_BOTH) {
         is_h_bond_donor_and_acceptor = true;
      }
   }

   if (is_h_bond_donor_and_acceptor == false) {
      if (resname_1 == "HOH") {
         if (hb_2 == HB_ACCEPTOR || hb_2 == HB_DONOR || hb_2 == HB_BOTH || hb_2 == HB_HYDROGEN) {
            is_h_bond_H_and_acceptor = true;
         }
      }
      if (resname_2 == "HOH") {
         if (hb_1 == HB_ACCEPTOR || hb_1 == HB_DONOR || hb_1 == HB_BOTH || hb_1 == HB_HYDROGEN) {
            is_h_bond_H_and_acceptor = true;
         }
      }
   }
}

std::string atom_overlaps_container_t::h_bond_info_t::format() const {
   return "h_bond_info_t";
}

} // namespace coot