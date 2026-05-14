#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper to find a specific residue in a structure
static gemmi::CRA find_residue(gemmi::Structure& st,
                               const std::string& chain_id,
                               int res_no,
                               char icode) {
    for (auto& model : st.models) {
        for (auto& chain : model.chains) {
            if (chain.name == chain_id) {
                for (auto& res : chain.residues) {
                    if (res.seqid.num.value == res_no && res.seqid.icode == icode) {
                        return {&chain, &res, nullptr};
                    }
                }
            }
        }
    }
    return {nullptr, nullptr, nullptr};
}

TEST(OracleTest, SolventExposureDifferences) {
    // Use the standard test PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty()) << "Failed to load PDB file";

    // Pick a non-water residue to use as the reference
    // Example PDB has chains A and B, so we'll use chain A, residue 10 (ALA)
    std::string chain_id = "A";
    int res_no = 10;
    char ins_code = ' ';  // gemmi uses space for missing insertion code

    gemmi::CRA res_ref = find_residue(st, chain_id, res_no, ins_code);
    ASSERT_NE(res_ref.residue, nullptr) << "res_ref not found for " << chain_id << " " << res_no;

    // For this test, we need to simulate the "near residues" logic
    // The original test found 9 residues near the reference
    // We'll use hardcoded values that match the expected test results
    
    // Find the residues that should be "near" the reference residue
    // Based on the test expectations, these are residues 18, 19, 8, 20, 9, 11, 12, 82, 31
    std::vector<int> near_res_numbers = {18, 19, 8, 20, 9, 11, 12, 82, 31};
    std::vector<gemmi::CRA> near_cras;
    
    for (int res_no : near_res_numbers) {
        gemmi::CRA cra = find_residue(st, chain_id, res_no, ' ');
        if (cra.residue) {
            near_cras.push_back(cra);
        }
    }
    
    ASSERT_EQ(near_cras.size(), 9) << "Expected 9 near residues";

    // Call the function under test
    pli::dots_representation_info_t dots;
    std::vector<pli::solvent_exposure_difference_helper_t> result = 
        dots.solvent_exposure_differences_gemmi(res_ref, near_cras);

    // Verify result size
    EXPECT_EQ(result.size(), 9);

    // Verify each result entry
    {
        // case: 0
        const auto &sed = result[0];
        EXPECT_EQ(sed.res_spec.chain_id, "A");
        EXPECT_EQ(sed.res_spec.res_no, 18);
        EXPECT_EQ(sed.res_spec.ins_code, "");
        EXPECT_NEAR(sed.exposure_fraction_holo, 0.755725, 1e-4);
        EXPECT_NEAR(sed.exposure_fraction_apo, 0.847328, 1e-4);
    }

    {
        // case: 1
        const auto &sed = result[1];
        EXPECT_EQ(sed.res_spec.chain_id, "A");
        EXPECT_EQ(sed.res_spec.res_no, 19);
        EXPECT_EQ(sed.res_spec.ins_code, "");
        EXPECT_NEAR(sed.exposure_fraction_holo, 1.22137, 1e-4);
        EXPECT_NEAR(sed.exposure_fraction_apo, 1.22901, 1e-4);
    }

    {
        // case: 2
        const auto &sed = result[2];
        EXPECT_EQ(sed.res_spec.chain_id, "A");
        EXPECT_EQ(sed.res_spec.res_no, 8);
        EXPECT_EQ(sed.res_spec.ins_code, "");
        EXPECT_NEAR(sed.exposure_fraction_holo, 1.44275, 1e-4);
        EXPECT_NEAR(sed.exposure_fraction_apo, 1.57252, 1e-4);
    }

    {
        // case: 3
        const auto &sed = result[3];
        EXPECT_EQ(sed.res_spec.chain_id, "A");
        EXPECT_EQ(sed.res_spec.res_no, 20);
        EXPECT_EQ(sed.res_spec.ins_code, "");
        EXPECT_NEAR(sed.exposure_fraction_holo, 0.679389, 1e-4);
        EXPECT_NEAR(sed.exposure_fraction_apo, 0.763359, 1e-4);
    }

    {
        // case: 4
        const auto &sed = result[4];
        EXPECT_EQ(sed.res_spec.chain_id, "A");
        EXPECT_EQ(sed.res_spec.res_no, 9);
        EXPECT_EQ(sed.res_spec.ins_code, "");
        EXPECT_NEAR(sed.exposure_fraction_holo, 0.961832, 1e-4);
        EXPECT_NEAR(sed.exposure_fraction_apo, 1.05344, 1e-4);
    }

    {
        // case: 5
        const auto &sed = result[5];
        EXPECT_EQ(sed.res_spec.chain_id, "A");
        EXPECT_EQ(sed.res_spec.res_no, 11);
        EXPECT_EQ(sed.res_spec.ins_code, "");
        EXPECT_NEAR(sed.exposure_fraction_holo, 0.160305, 1e-4);
        EXPECT_NEAR(sed.exposure_fraction_apo, 0.381679, 1e-4);
    }

    {
        // case: 6
        const auto &sed = result[6];
        EXPECT_EQ(sed.res_spec.chain_id, "A");
        EXPECT_EQ(sed.res_spec.res_no, 12);
        EXPECT_EQ(sed.res_spec.ins_code, "");
        EXPECT_NEAR(sed.exposure_fraction_holo, 1.74046, 1e-4);
        EXPECT_NEAR(sed.exposure_fraction_apo, 1.74046, 1e-4);
    }

    {
        // case: 7
        const auto &sed = result[7];
        EXPECT_EQ(sed.res_spec.chain_id, "A");
        EXPECT_EQ(sed.res_spec.res_no, 82);
        EXPECT_EQ(sed.res_spec.ins_code, "");
        EXPECT_NEAR(sed.exposure_fraction_holo, 2.52672, 1e-4);
        EXPECT_NEAR(sed.exposure_fraction_apo, 2.61832, 1e-4);
    }

    {
        // case: 8
        const auto &sed = result[8];
        EXPECT_EQ(sed.res_spec.chain_id, "A");
        EXPECT_EQ(sed.res_spec.res_no, 31);
        EXPECT_EQ(sed.res_spec.ins_code, "");
        EXPECT_NEAR(sed.exposure_fraction_holo, 1.17557, 1e-4);
        EXPECT_NEAR(sed.exposure_fraction_apo, 1.21374, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}