#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_HA_unit_vector) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    ASSERT_FALSE(st.models.empty());
    const auto& model = st.models[0];
    
    // Collect residues
    std::vector<const gemmi::Residue*> residues;
    for (const auto& chain : model.chains)
        for (const auto& res : chain.residues)
            residues.push_back(&res);
    
    ASSERT_GT(residues.size(), 0) << "Should have found residues";
    
    // Test case 1: residue 6 (has CA,C,N,CB) - should succeed
    {
        const gemmi::Residue* res = residues[5];  // 0-indexed, residue 6
        auto result = coot::molecule_t_get_HA_unit_vector_gemmi(*res);
        EXPECT_TRUE(result.first) << "Residue 6 should have all required atoms";
        if (result.first) {
            EXPECT_NEAR(result.second.x, 0.177931, 3e-6);
            EXPECT_NEAR(result.second.y, -0.635763, 3e-6);
            EXPECT_NEAR(result.second.z, -0.751097, 3e-6);
        }
    }
    
    // Test case 2: residue 11 (GLY, has CA,C,N but missing CB) - should succeed
    {
        const gemmi::Residue* res = residues[10];  // 0-indexed, residue 11
        auto result = coot::molecule_t_get_HA_unit_vector_gemmi(*res);
        EXPECT_TRUE(result.first) << "GLY 11 should have required atoms (CA,C,N) for get_HA_unit_vector";
        if (result.first) {
            EXPECT_NEAR(result.second.x, 0.747214, 3e-6);
            EXPECT_NEAR(result.second.y, 0.315477, 3e-6);
            EXPECT_NEAR(result.second.z, -0.584932, 3e-6);
        }
    }
    
    // Test case 3: last residue (HOH 1 in chain B, only has O atom) - should fail
    {
        const gemmi::Residue* res = residues.back();
        auto result = coot::molecule_t_get_HA_unit_vector_gemmi(*res);
        EXPECT_FALSE(result.first) << "Water residue should not have required atoms";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}