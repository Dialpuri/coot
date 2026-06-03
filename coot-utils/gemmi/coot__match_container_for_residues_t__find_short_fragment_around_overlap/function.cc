#include "function.hh"

namespace coot {

delete_a_chain_t match_container_for_residues_t::find_short_fragment_around_overlap_gemmi(
    const std::vector<gemmi::Atom*>& atom_selection_1,
    const std::vector<gemmi::Atom*>& atom_selection_2) const {

   // these fragments are *not* mergeable if the same fragment is short on *both* sides.
   // Consider the case of a fragment 1-20 and a fragment that matches 5-12. We don't want
   // to merge that - just delete the short fragment.

   bool is_mergeable = true;

   bool is_first_selection = true;
   bool is_upstream = true;

   // how many atoms below the matched atoms?
   int atom_sel_1_n_above = 0;
   int atom_sel_1_n_below = 0;
   int atom_sel_2_n_above = 0;
   int atom_sel_2_n_below = 0;

   bool hit_matches = false;
   for (int iat = 0; iat < static_cast<int>(atom_selection_1.size()); iat++) {
      gemmi::Atom *at = atom_selection_1[iat];
      if (!hit_matches) {
         for (unsigned int ip = 0; ip < atom_pairs.size(); ip++) {
            if (atom_pairs[ip].first == at) {
               hit_matches = true;
               break;
            }
         }
      }
      if (!hit_matches)
         atom_sel_1_n_above++;
      else
         atom_sel_1_n_below++;
   }
   hit_matches = false;
   for (int iat = 0; iat < static_cast<int>(atom_selection_2.size()); iat++) {
      gemmi::Atom *at = atom_selection_2[iat];
      if (!hit_matches) {
         for (unsigned int ip = 0; ip < atom_pairs.size(); ip++) {
            if (atom_pairs[ip].second == at) {
               hit_matches = true;
               break;
            }
         }
      }
      if (!hit_matches)
         atom_sel_2_n_above++;
      else
         atom_sel_2_n_below++;
   }

   // idx == 1: is_in_first = true;  is_upstream = true
   // idx == 2: is_in_first = true;  is_upstream = false
   // idx == 3: is_in_first = false; is_upstream = true
   // idx == 4: is_in_first = false; is_upstream = false

   delete_a_chain_t::delete_a_chain_type_t dt = delete_a_chain_t::NONE;

   // these fragments should *not* mergeable if the same fragment is short on *both* sides.
   // (or as tested here) the same fragment is long on both sides
   if (atom_sel_1_n_below >= atom_sel_2_n_below) {
      if (atom_sel_1_n_above >= atom_sel_2_n_above) {
         is_mergeable = false;
         dt = delete_a_chain_t::DELETE_SECOND_CHAIN;
      }
   }
   if (atom_sel_2_n_below >= atom_sel_1_n_below) {
      if (atom_sel_2_n_above >= atom_sel_1_n_above) {
         is_mergeable = false;
         dt = delete_a_chain_t::DELETE_FIRST_CHAIN;
      }
   }

   int idx = 1;
   if (atom_sel_1_n_below <= atom_sel_1_n_above)
      if (atom_sel_1_n_below <= atom_sel_2_n_above)
         if (atom_sel_1_n_below <= atom_sel_2_n_below)
            idx = 2;
   if (atom_sel_2_n_above <= atom_sel_1_n_above)
      if (atom_sel_2_n_above <= atom_sel_1_n_below)
         if (atom_sel_2_n_above <= atom_sel_2_n_below)
            idx = 3;
   if (atom_sel_2_n_below <= atom_sel_1_n_above)
      if (atom_sel_2_n_below <= atom_sel_1_n_below)
         if (atom_sel_2_n_below <= atom_sel_2_n_above)
            idx = 4;

   if (idx == 2) { is_upstream = false; }
   if (idx == 3) { is_first_selection = false; }
   if (idx == 4) { is_first_selection = false; is_upstream = false; }

   delete_a_chain_t dac(is_mergeable, is_first_selection, is_upstream);
   dac.delete_type = dt;
   return dac;
}

} // namespace coot