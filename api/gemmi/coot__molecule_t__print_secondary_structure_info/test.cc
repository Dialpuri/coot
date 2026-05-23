#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, print_secondary_structure_info) {
    // case: example.pdb with no secondary structure
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        EXPECT_FALSE(st.models.empty());

        // Placeholder: just ensure the function compiles and runs without throwing
        coot::molecule_t mol;
        EXPECT_NO_THROW(mol.print_secondary_structure_info_gemmi());
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}