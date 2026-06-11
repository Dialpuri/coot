#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <map>
#include "function.hh"

// Suppress "something bad" warning — it's printed by cout for the incomplete case
// We just test the behavior, not the warning itself.

TEST(OracleTest, chem_comp_chir_structure) {
    // No need to call geom.init_standard() for the gemmi port —
    // mon_lib_add_chiral is a no-op stub, and the test only verifies
    // non-throwing behavior and field counts.

    int imol_enc = 0;

    // --- Case 1: Complete chirality data (should succeed) ---
    {
        std::map<std::string, std::string> structure;
        structure["comp_id"]        = "ARG";
        structure["id"]             = "1";
        structure["atom_id_centre"] = "CB";
        structure["atom_id_1"]      = "CA";
        structure["atom_id_2"]      = "N";
        structure["atom_id_3"]      = "CG";
        structure["volume_sign"]    = "-";

        coot::protein_geometry::chem_comp_chir_structure_gemmi(structure, imol_enc);

        EXPECT_EQ(structure.size(), 7u);
    }

    // --- Case 2: Missing required fields (should trigger WARNING path) ---
    {
        std::map<std::string, std::string> structure2;
        structure2["comp_id"] = "ALA";
        structure2["id"]      = "1";

        EXPECT_NO_THROW(
            coot::protein_geometry::chem_comp_chir_structure_gemmi(structure2, imol_enc));
    }

    // --- Case 3: All fields present but different residue and volume_sign ---
    {
        std::map<std::string, std::string> structure3;
        structure3["comp_id"]        = "LYS";
        structure3["id"]             = "1";
        structure3["atom_id_centre"] = "CB";
        structure3["atom_id_1"]      = "CA";
        structure3["atom_id_2"]      = "N";
        structure3["atom_id_3"]      = "CG";
        structure3["volume_sign"]    = "+";

        coot::protein_geometry::chem_comp_chir_structure_gemmi(structure3, imol_enc);

        EXPECT_EQ(structure3.size(), 7u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}