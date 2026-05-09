#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, FillPartialResidues) {
   std::string pdb_path = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb";
   gemmi::Structure st = gemmi::read_pdb_file(pdb_path);
   
   int status = coot::fill_partial_residues_gemmi(st);
   EXPECT_EQ(status, 0);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}