#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, copy_chain) {
    // Load the test PDB file
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_EQ(st.models.size(), 1) << "Expected exactly one model in PDB file";
    
    // Test case 1: copy chain A
    {
        gemmi::Model& model = st.models[0];
        gemmi::Chain* chain_p = nullptr;
        for (auto& c : model.chains) {
            if (c.name == "A") {
                chain_p = &c;
                break;
            }
        }
        
        ASSERT_NE(chain_p, nullptr) << "Chain A not found";
        EXPECT_EQ(chain_p->name, "A");
        EXPECT_EQ(chain_p->residues.size(), 267);
        
        std::pair<gemmi::Chain*, gemmi::Structure> result = coot::util::copy_chain_gemmi(*chain_p);
        gemmi::Chain* new_chain = result.first;
        gemmi::Structure& new_st = result.second;
        
        EXPECT_TRUE(new_chain != nullptr);
        EXPECT_EQ(new_st.models.size(), 1);
        
        if (new_chain && !new_st.models.empty()) {
            gemmi::Model& new_model = new_st.models[0];
            EXPECT_EQ(new_chain->name, "A");
            EXPECT_EQ(new_chain->residues.size(), 267);
            EXPECT_EQ(new_model.chains.size(), 1);
        }
    }
    
    // Test case 2: copy chain B
    {
        gemmi::Model& model = st.models[0];
        gemmi::Chain* chain_p = nullptr;
        for (auto& c : model.chains) {
            if (c.name == "B") {
                chain_p = &c;
                break;
            }
        }
        
        ASSERT_NE(chain_p, nullptr) << "Chain B not found";
        EXPECT_EQ(chain_p->name, "B");
        EXPECT_EQ(chain_p->residues.size(), 1);
        
        std::pair<gemmi::Chain*, gemmi::Structure> result = coot::util::copy_chain_gemmi(*chain_p);
        gemmi::Chain* new_chain = result.first;
        gemmi::Structure& new_st = result.second;
        
        EXPECT_TRUE(new_chain != nullptr);
        EXPECT_EQ(new_st.models.size(), 1);
        
        if (new_chain && !new_st.models.empty()) {
            gemmi::Model& new_model = new_st.models[0];
            EXPECT_EQ(new_chain->name, "B");
            EXPECT_EQ(new_chain->residues.size(), 1);
            EXPECT_EQ(new_model.chains.size(), 1);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}