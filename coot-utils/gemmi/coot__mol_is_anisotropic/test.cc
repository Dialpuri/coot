#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, mol_is_anisotropic) {
    // Case 1: Test with example.pdb (likely isotropic)
    { // case: example.pdb
        gemmi::Structure st1 = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        bool is_aniso = coot::mol_is_anisotropic_gemmi(st1);
        EXPECT_FALSE(is_aniso);
    }
    
    // Case 2: Load example-hydrogen.pdb (also likely isotropic)
    { // case: example-hydrogen.pdb
        gemmi::Structure st2 = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
        bool is_aniso = coot::mol_is_anisotropic_gemmi(st2);
        EXPECT_FALSE(is_aniso);
    }
    
    // Case 3: Test with empty structure (graceful handling)
    { // case: empty structure
        gemmi::Structure st_empty;
        bool is_aniso_empty = coot::mol_is_anisotropic_gemmi(st_empty);
        EXPECT_FALSE(is_aniso_empty);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}