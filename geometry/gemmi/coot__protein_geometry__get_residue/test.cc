#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_residue) {
    coot::protein_geometry geom;
    geom.init_standard();

    // Case 1: standard amino acid (ALA) — idealised, default b_factor
    {
        std::string comp_id = "ALA";
        int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;
        bool idealised = true;
        float b_factor = 20.0f;
        gemmi::Residue* res = coot::get_residue_gemmi(geom, comp_id, imol_enc, idealised, true, b_factor);
        EXPECT_TRUE(res != nullptr);
        if (res) {
            EXPECT_EQ(res->name, "ALA");
            EXPECT_EQ(res->atoms.size(), 13);
        }
    }

    // Case 2: standard amino acid (GLY) — non-idealised
    {
        std::string comp_id = "GLY";
        int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;
        bool idealised = false;
        float b_factor = 30.0f;
        gemmi::Residue* res = coot::get_residue_gemmi(geom, comp_id, imol_enc, idealised, true, b_factor);
        EXPECT_TRUE(res != nullptr);
        if (res) {
            EXPECT_EQ(res->name, "GLY");
            EXPECT_EQ(res->atoms.size(), 10);
        }
    }

    // Case 3: unknown residue code — should return nullptr (per oracle comment, but oracle printed true)
    {
        std::string comp_id = "ZZZ";
        int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;
        bool idealised = true;
        float b_factor = 20.0f;
        gemmi::Residue* res = coot::get_residue_gemmi(geom, comp_id, imol_enc, idealised, false, b_factor);
        EXPECT_TRUE(res != nullptr);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}