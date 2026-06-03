#include "function.hh"
#include <vector>
#include <set>
#include "geometry/residue-and-atom-specs.hh"

namespace coot {

void restraints_container_t::fill_old_to_new_index_vector() {
}

bool restraints_container_t::try_add_using_old_atom_indices(const extra_restraints_t::extra_angle_restraint_t& ebr) {
   return false;
}

bool restraints_container_t::fixed_check(int idx) {
   return false;
}

std::vector<bool> restraints_container_t::make_fixed_flags(int idx_1, int idx_2, int idx_3) {
   return std::vector<bool>(3, false);
}

void restraints_container_t::add_user_defined_angle_restraint(
   int type, int idx_1, int idx_2, int idx_3,
   const std::vector<bool>& fixed_flags,
   double angle, double esd, double dummy) {
   std::vector<int> restraint;
   restraint.push_back(type);
   restraint.push_back(idx_1);
   restraint.push_back(idx_2);
   restraint.push_back(idx_3);
   restraints_vec.push_back(restraint);
}

void add_extra_angle_restraints_gemmi(restraints_container_t& rc,
                                     const extra_restraints_t& extra_restraints) {
   rc.fill_old_to_new_index_vector();

   for (unsigned int i = 0; i < extra_restraints.angle_restraints.size(); i++) {
      const extra_restraints_t::extra_angle_restraint_t& ebr = extra_restraints.angle_restraints[i];

      bool done = rc.try_add_using_old_atom_indices(ebr);
      if (done) continue;

      if (rc.from_residue_vector) {
      } else {
      }
   }
}

} // namespace coot