#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

namespace coot {

TEST(AddExtraTorsionRestraints, PortedFunction) {
   restraints_container_t rc;
   rc.from_residue_vector = false;
   rc.udd_atom_index_handle = 0;
   rc.mol = nullptr;
   
   extra_restraints_t extra_restraints;
   add_extra_torsion_restraints_gemmi(rc, extra_restraints);
}

TEST(AddExtraTorsionRestraints, WithRestraints) {
   restraints_container_t rc;
   rc.from_residue_vector = false;
   rc.udd_atom_index_handle = 0;
   rc.mol = nullptr;
   
   extra_restraints_t extra_restraints;
   atom_spec_t spec_1("A", 10, " ", "CA", "");
   atom_spec_t spec_2("A", 10, " ", "C",  "");
   atom_spec_t spec_3("A", 11, " ", "N",  "");
   atom_spec_t spec_4("A", 11, " ", "CA", "");
   extra_restraints.torsion_restraints.push_back(
      extra_restraints_t::extra_torsion_restraint_t(
         spec_1, spec_2, spec_3, spec_4, -120.0, 20.0, 1)
   );
   add_extra_torsion_restraints_gemmi(rc, extra_restraints);
}

} // namespace coot

int main(int argc, char** argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}