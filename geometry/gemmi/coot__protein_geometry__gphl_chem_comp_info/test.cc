#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, gphl_chem_comp_info) {
    coot::protein_geometry geom;
    geom.init_standard();

    // Case 1: Struct with comp_id = "ALA" (present in standard dict)
    {
        std::map<std::string, std::string> structure_data;
        structure_data["comp_id"] = "ALA";
        structure_data["arguments"] = "test_args";
        structure_data["run_date"] = "2024-01-01";

        int imol_enc = 0;

        // Record state before: find ALA's gphl info size
        int gphl_size_before = -1;
        int idx = geom.get_monomer_restraints_index("ALA", imol_enc, true);
        if (idx >= 0) {
            gphl_size_before = static_cast<int>(geom[idx].second.gphl_chem_comp_info.info.size());
        }

        // Call the ported function
        coot::gphl_chem_comp_info_t result = coot::gphl_chem_comp_info_gemmi(geom, structure_data, imol_enc);

        // Record state after: use the return value size
        int gphl_size_after = static_cast<int>(result.info.size());

        EXPECT_EQ(gphl_size_before, 0);
        EXPECT_EQ(gphl_size_after, 3);
    }

    // Case 2: Struct without comp_id — should not add anything
    {
        std::map<std::string, std::string> structure_data;
        structure_data["arguments"] = "no_comp_id_here";

        int imol_enc = 0;

        int dict_size_before = static_cast<int>(geom.size());

        coot::gphl_chem_comp_info_gemmi(geom, structure_data, imol_enc);

        int dict_size_after = static_cast<int>(geom.size());

        EXPECT_EQ(dict_size_before, 41);
        EXPECT_EQ(dict_size_after, 41);
    }

    // Case 3: Struct with comp_id for unknown residue — should not match dict_res_restraints
    {
        std::map<std::string, std::string> structure_data;
        structure_data["comp_id"] = "ZZZ";
        structure_data["arguments"] = "unknown_ligand";

        int imol_enc = 0;

        int dict_size_before = static_cast<int>(geom.size());

        coot::gphl_chem_comp_info_gemmi(geom, structure_data, imol_enc);

        int dict_size_after = static_cast<int>(geom.size());

        EXPECT_EQ(dict_size_before, 41);
        EXPECT_EQ(dict_size_after, 41);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}