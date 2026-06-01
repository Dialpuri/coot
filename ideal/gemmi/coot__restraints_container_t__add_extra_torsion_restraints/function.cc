#include "function.hh"
#include <vector>
#include <set>
#include "geometry/residue-and-atom-specs.hh"

namespace coot {

void restraints_container_t::fill_old_to_new_index_vector() {
}

bool restraints_container_t::try_add_using_old_atom_indices(const extra_restraints_t::extra_torsion_restraint_t& etr) {
   return false;
}

bool restraints_container_t::fixed_check(int idx) {
   return false;
}

std::vector<bool> restraints_container_t::make_fixed_flags(int idx_1, int idx_2, int idx_3, int idx_4) {
   return std::vector<bool>(4, false);
}

void restraints_container_t::add_torsion(int type, int idx_1, int idx_2, int idx_3, int idx_4,
                                          const std::vector<bool>& fixed_flags,
                                          double angle, double esd, double dummy, int period) {
   std::vector<int> restraint;
   restraint.push_back(type);
   restraint.push_back(idx_1);
   restraint.push_back(idx_2);
   restraint.push_back(idx_3);
   restraint.push_back(idx_4);
   restraints_vec.push_back(restraint);
}

void add_extra_torsion_restraints_gemmi(restraints_container_t& rc,
                                         const extra_restraints_t& extra_restraints) {
   rc.fill_old_to_new_index_vector();

   for (unsigned int i = 0; i < extra_restraints.torsion_restraints.size(); i++) {
      const extra_restraints_t::extra_torsion_restraint_t& etr = extra_restraints.torsion_restraints[i];

      bool done = rc.try_add_using_old_atom_indices(etr);
      if (done) continue;

      if (rc.from_residue_vector) {
      } else {
      }
   }
}

} // namespace coot