#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

namespace coot {

TEST(AddExtraBondRestraints, PortedFunction) {
   restraints_container_t rc;
   rc.from_residue_vector = false;
   rc.udd_atom_index_handle = 0;
   
   extra_restraints_t extra_restraints;
   add_extra_bond_restraints_gemmi(rc, extra_restraints);
}

TEST(AddExtraBondRestraints, WithRestraints) {
   restraints_container_t rc;
   rc.from_residue_vector = false;
   rc.udd_atom_index_handle = 0;
   
   extra_restraints_t extra_restraints;
   atom_spec_t spec_1("A", 10, " ", "CA", "");
   atom_spec_t spec_2("A", 11, " ", "CA", "");
   extra_restraints.bond_restraints.push_back(
      extra_restraints_t::extra_bond_restraint_t(spec_1, spec_2, 3.8, 0.1)
   );
   add_extra_bond_restraints_gemmi(rc, extra_restraints);
}

} // namespace coot

int main(int argc, char** argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}