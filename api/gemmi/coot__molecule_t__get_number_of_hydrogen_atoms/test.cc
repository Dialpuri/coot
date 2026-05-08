#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_number_of_hydrogen_atoms_gemmi) {
  std::string pdb_path = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb";
  gemmi::Structure st = gemmi::read_pdb_file(pdb_path);
  coot::molecule_t mol;
  int n_h = mol.get_number_of_hydrogen_atoms_gemmi(st);
  EXPECT_EQ(n_h, 2173);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}