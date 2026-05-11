#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, ClosestApproach) {
    // Load the structure
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty()) << "Failed to load PDB file";
    
    gemmi::Model& model = st.models[0];
    
    // Helper lambda to get residue by chain ID and sequence number
    auto get_residue = [](gemmi::Model& model, const std::string& chain_name, int seqnum) -> gemmi::Residue* {
        for (auto& chain : model.chains) {
            if (chain.name == chain_name) {
                for (auto& res : chain.residues) {
                    if (res.seqid.num.value == seqnum) {
                        return &res;
                    }
                }
            }
        }
        return nullptr;
    };
    
    // case 1: two residues that exist and are close in sequence (spatially close)
    {
        gemmi::Residue* r1 = get_residue(model, "A", 10);
        gemmi::Residue* r2 = get_residue(model, "A", 12);
        
        ASSERT_NE(r1, nullptr) << "get_residue failed for r1 //A/10";
        ASSERT_NE(r2, nullptr) << "get_residue failed for r2 //A/12";
        
        std::pair<bool,float> result = coot::closest_approach_gemmi(*r1, *r2);
        
        EXPECT_TRUE(result.first);
        EXPECT_NEAR(result.second, 4.80432, 1e-4);
    }
    
    // case 2: two distant residues (sequence-wise and spatially)
    {
        gemmi::Residue* r1 = get_residue(model, "A", 10);
        gemmi::Residue* r2 = get_residue(model, "A", 50);
        
        ASSERT_NE(r1, nullptr) << "get_residue failed for r1 //A/10";
        ASSERT_NE(r2, nullptr) << "get_residue failed for r2 //A/50";
        
        std::pair<bool,float> result = coot::closest_approach_gemmi(*r1, *r2);
        
        EXPECT_TRUE(result.first);
        EXPECT_NEAR(result.second, 23.6919, 1e-4);
    }
    
    // case 3: same residue (edge case)
    {
        gemmi::Residue* r1 = get_residue(model, "A", 10);
        gemmi::Residue* r2 = get_residue(model, "A", 10);
        
        ASSERT_NE(r1, nullptr) << "get_residue failed for r1 //A/10";
        ASSERT_NE(r2, nullptr) << "get_residue failed for r2 //A/10";
        
        std::pair<bool,float> result = coot::closest_approach_gemmi(*r1, *r2);
        
        EXPECT_TRUE(result.first);
        EXPECT_NEAR(result.second, 0.0, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}