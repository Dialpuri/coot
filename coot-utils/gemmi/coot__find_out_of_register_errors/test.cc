#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, FindOutOfRegisterErrors) {
    // Read the PDB files
    gemmi::Structure ref_st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Structure mut_st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // The function prints output but returns void
    // We verify the function executes without throwing exceptions
    EXPECT_NO_THROW(coot::find_out_of_register_errors_gemmi(&mut_st, &ref_st));
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}