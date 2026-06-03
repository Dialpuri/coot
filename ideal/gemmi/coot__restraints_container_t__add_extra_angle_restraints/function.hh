#pragma once
#include <vector>
#include <string>
#include <set>
#include <utility>
#include "ideal/extra-restraints.hh"
#include "geometry/residue-and-atom-specs.hh"

namespace coot {

struct restraints_container_t {
   std::vector<std::pair<bool, void*>> residues_vec;
   bool from_residue_vector;
   std::vector<std::vector<int>> restraints_vec;
   std::vector<std::set<int>> bonded_atom_indices;
   int udd_atom_index_handle;
   void* mol;

   void fill_old_to_new_index_vector();
   bool try_add_using_old_atom_indices(const extra_restraints_t::extra_angle_restraint_t& ebr);
   bool fixed_check(int idx);
   std::vector<bool> make_fixed_flags(int idx_1, int idx_2, int idx_3);
   void add_user_defined_angle_restraint(
      int type, int idx_1, int idx_2, int idx_3,
      const std::vector<bool>& fixed_flags,
      double angle, double esd, double dummy);
};

void add_extra_angle_restraints_gemmi(restraints_container_t& rc,
                                     const extra_restraints_t& extra_restraints);

} // namespace coot