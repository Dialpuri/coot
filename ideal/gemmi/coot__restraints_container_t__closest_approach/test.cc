#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper function to get residue by chain ID and residue number from a gemmi Structure
static gemmi::Residue*
get_residue_by_seqid(gemmi::Structure& st, const std::string& chain_name, int seqnum) {
    for (auto& model : st.models)
        for (auto& chain : model.chains)
            if (chain.name == chain_name)
                for (auto& res : chain.residues)
                    if (res.seqid.num.value == seqnum)
                        return &res;
    return nullptr;
}

TEST(OracleTest, ClosestApproachTest) {
    // case: normal case - two different residues
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        
        gemmi::Residue *r1 = get_residue_by_seqid(st, "A", 10);
        gemmi::Residue *r2 = get_residue_by_seqid(st, "A", 20);
        
        ASSERT_NE(r1, nullptr) << "Could not find residue A/10";
        ASSERT_NE(r2, nullptr) << "Could not find residue A/20";
        
        // Test closest_approach
        std::pair<bool, float> result = coot::closest_approach_gemmi(*r1, *r2);
        
        EXPECT_TRUE(result.first) << "closest_approach should find a result";
        EXPECT_NEAR(result.second, 3.65605, 1e-4) << "Distance should match expected value";
    }

    // case: edge case - same residue
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        
        gemmi::Residue *r1 = get_residue_by_seqid(st, "A", 10);
        ASSERT_NE(r1, nullptr) << "Could not find residue A/10";
        
        // Test closest_approach with same residue
        std::pair<bool, float> result_same = coot::closest_approach_gemmi(*r1, *r1);
        
        EXPECT_TRUE(result_same.first) << "closest_approach should find a result for same residue";
        EXPECT_NEAR(result_same.second, 0.0, 1e-4) << "Distance should be 0 for same residue";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}