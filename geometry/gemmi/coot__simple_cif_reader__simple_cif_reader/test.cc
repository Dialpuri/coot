#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, simple_cif_reader) {
    // Case 1: Non-existent file — should trigger warning and empty vectors
    {
        std::string cif_dict = "/nonexistent/chem_comp.cif";
        coot::simple_cif_reader reader;
        reader.simple_cif_reader_gemmi(cif_dict);
        EXPECT_EQ(reader.names.size(), 0u);
        EXPECT_EQ(reader.three_letter_codes.size(), 0u);
    }

    // Case 2: Try a path that might exist (the test cif file)
    {
        std::string cif_dict = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-ligand.cif";
        coot::simple_cif_reader reader;
        reader.simple_cif_reader_gemmi(cif_dict);
        EXPECT_EQ(reader.names.size(), 1u);
        EXPECT_EQ(reader.three_letter_codes.size(), 1u);
        if (!reader.three_letter_codes.empty()) {
            EXPECT_EQ(reader.three_letter_codes[0], "LZA");
        }
        // Test has_restraints_for with a common residue
        bool has_val = reader.has_restraints_for("ALA");
        EXPECT_EQ(has_val, false);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}