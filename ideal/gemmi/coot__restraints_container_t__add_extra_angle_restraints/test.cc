#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

namespace coot {

TEST(AddExtraAngleRestraints, BasicFunctionCall) {
   restraints_container_t rc;
   rc.from_residue_vector = false;
   rc.udd_atom_index_handle = 0;
   
   extra_restraints_t extra_restraints;
   add_extra_angle_restraints_gemmi(rc, extra_restraints);
}

TEST(AddExtraAngleRestraints, WithSingleAngleRestraint) {
   restraints_container_t rc;
   rc.from_residue_vector = false;
   rc.udd_atom_index_handle = 0;
   
   extra_restraints_t extra_restraints;
   atom_spec_t spec_1("A", 5, " ", "N", "");
   atom_spec_t spec_2("A", 5, " ", "CA", "");
   atom_spec_t spec_3("A", 5, " ", "C", "");
   extra_restraints.angle_restraints.push_back(
      extra_restraints_t::extra_angle_restraint_t(spec_1, spec_2, spec_3, 120.0, 5.0)
   );
   add_extra_angle_restraints_gemmi(rc, extra_restraints);
}

TEST(AddExtraAngleRestraints, WithInvalidAtoms) {
   restraints_container_t rc;
   rc.from_residue_vector = false;
   rc.udd_atom_index_handle = 0;
   
   extra_restraints_t extra_restraints;
   atom_spec_t spec_1("A", 99999, " ", "N", "");
   atom_spec_t spec_2("A", 99999, " ", "CA", "");
   atom_spec_t spec_3("A", 99999, " ", "C", "");
   extra_restraints.angle_restraints.push_back(
      extra_restraints_t::extra_angle_restraint_t(spec_1, spec_2, spec_3, 120.0, 5.0)
   );
   add_extra_angle_restraints_gemmi(rc, extra_restraints);
}

TEST(AddExtraAngleRestraints, WithMultipleAngleRestraints) {
   restraints_container_t rc;
   rc.from_residue_vector = false;
   rc.udd_atom_index_handle = 0;
   
   extra_restraints_t extra_restraints;
   
   // First angle: CA-C-O in residue 5
   atom_spec_t spec_1_1("A", 5, " ", "CA", "");
   atom_spec_t spec_1_2("A", 5, " ", "C", "");
   atom_spec_t spec_1_3("A", 5, " ", "O", "");
   extra_restraints.angle_restraints.push_back(
      extra_restraints_t::extra_angle_restraint_t(spec_1_1, spec_1_2, spec_1_3, 120.0, 5.0)
   );
   
   // Second angle: C-N-CA in residues 5-6
   atom_spec_t spec_2_1("A", 5, " ", "C", "");
   atom_spec_t spec_2_2("A", 6, " ", "N", "");
   atom_spec_t spec_2_3("A", 6, " ", "CA", "");
   extra_restraints.angle_restraints.push_back(
      extra_restraints_t::extra_angle_restraint_t(spec_2_1, spec_2_2, spec_2_3, 125.0, 5.0)
   );
   
   add_extra_angle_restraints_gemmi(rc, extra_restraints);
}

} // namespace coot

int main(int argc, char** argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}