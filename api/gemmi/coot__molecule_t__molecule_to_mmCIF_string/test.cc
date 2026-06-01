#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, molecule_to_mmCIF_string) {
    // Case 1: valid protein PDB
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        EXPECT_FALSE(st.models.empty());

        std::string cif = coot::molecule_to_mmCIF_string_gemmi(st);
        EXPECT_FALSE(cif.empty());

        // Check preview: should contain mmCIF markers
        EXPECT_NE(cif.find("data_"), std::string::npos);
        EXPECT_NE(cif.find("_entry"), std::string::npos);
        EXPECT_NE(cif.find("_struct"), std::string::npos);
    }

    // Case 2: hydrogenated PDB (should also produce CIF)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
        EXPECT_FALSE(st.models.empty());

        std::string cif = coot::molecule_to_mmCIF_string_gemmi(st);
        EXPECT_FALSE(cif.empty());
    }

    // Case 3: invalid molecule (edge case) - empty structure
    {
        gemmi::Structure st;
        std::string cif = coot::molecule_to_mmCIF_string_gemmi(st);
        EXPECT_TRUE(cif.empty());
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}