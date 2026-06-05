#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, inverted_chiral_volumes) {
    // Initialize protein geometry (required pattern)
    coot::protein_geometry geom;
    geom.init_standard();

    int imol = 0;
    int cif_dictionary_read_number = 0;

    // Case 1: example.pdb
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        ASSERT_GT(st.models.size(), 0u);

        std::pair<std::vector<std::string>, std::vector<coot::atom_spec_t>> result =
            coot::inverted_chiral_volumes_gemmi(imol, st, &geom, cif_dictionary_read_number);

        // Oracle output: unknown_residue_types: 0, inverted_chirals_count: 0
        EXPECT_EQ(result.first.size(), 0u);
        EXPECT_EQ(result.second.size(), 0u);
    }

    // Case 2: example-hydrogen.pdb
    {
        gemmi::Structure st2 = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

        ASSERT_GT(st2.models.size(), 0u);

        std::pair<std::vector<std::string>, std::vector<coot::atom_spec_t>> result2 =
            coot::inverted_chiral_volumes_gemmi(imol, st2, &geom, cif_dictionary_read_number);

        // Oracle output: unknown_residue_types: 0, inverted_chirals_count: 0
        EXPECT_EQ(result2.first.size(), 0u);
        EXPECT_EQ(result2.second.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}