#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, chem_comp_component) {
    coot::protein_geometry geom;
    geom.init_standard();

    // Case 1: Complete _chem_comp structure with all required fields
    {
        std::map<std::string, std::string> fields;
        fields["id"] = "LEU";
        fields["three_letter_code"] = "LEU";
        fields["name"] = "LEUCINE";
        fields["type"] = "L-Peptide";
        fields["number_atoms_all"] = "20";
        fields["number_atoms_nh"] = "14";

        int imol_enc = 0;
        EXPECT_EQ(fields.size(), 6);
        EXPECT_EQ(imol_enc, 0);

        std::string result = coot::chem_comp_component_gemmi(fields, imol_enc);
        EXPECT_EQ(result, "LEU");
        EXPECT_FALSE(result.empty());
    }

    // Case 2: Missing required field (no three_letter_code)
    {
        std::map<std::string, std::string> fields;
        fields["id"] = "ALA";
        fields["name"] = "ALANINE";
        // Missing: three_letter_code

        int imol_enc = 0;
        EXPECT_EQ(fields.size(), 2);
        EXPECT_EQ(imol_enc, 0);

        std::string result = coot::chem_comp_component_gemmi(fields, imol_enc);
        EXPECT_EQ(result, "ALA");
        EXPECT_FALSE(result.empty());
    }

    // Case 3: Missing 'id' field entirely
    {
        std::map<std::string, std::string> fields;
        fields["three_letter_code"] = "GLY";
        fields["name"] = "GLYCINE";
        // Missing: id

        int imol_enc = 0;
        EXPECT_EQ(fields.size(), 2);
        EXPECT_EQ(imol_enc, 0);

        std::string result = coot::chem_comp_component_gemmi(fields, imol_enc);
        EXPECT_EQ(result, "");
        EXPECT_TRUE(result.empty());
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}