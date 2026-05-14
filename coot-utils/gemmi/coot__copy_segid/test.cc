#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <cstring>
#include <string>
#include "function.hh"

TEST(OracleTest, copy_segid) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Find residues by chain ID and sequence number
    gemmi::Residue* provider = nullptr;
    gemmi::Residue* receiver = nullptr;
    
    for (auto& model : st.models) {
        for (auto& chain : model.chains) {
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == 10 && chain.name == "A") {
                    provider = &res;
                }
                if (res.seqid.num.value == 11 && chain.name == "A") {
                    receiver = &res;
                }
            }
        }
    }
    
    ASSERT_NE(provider, nullptr) << "Could not find provider residue (A/10)";
    ASSERT_NE(receiver, nullptr) << "Could not find receiver residue (A/11)";
    
    // Case 1: copy segid from a provider residue to a receiver residue
    {
        // Modify provider segid to something non-empty
        std::string new_segid = "PROV";
        provider->segment = new_segid;
        
        // Capture receiver segid before
        std::string before_segid = receiver->segment;
        
        // Call the function
        bool result = coot::copy_segid_gemmi(*provider, *receiver);
        EXPECT_TRUE(result);
        
        // Verify receiver segid after
        EXPECT_EQ(receiver->segment, "PROV");
    }
    
    // Case 2: test with empty segid provider (edge case)
    {
        gemmi::Residue* provider2 = nullptr;
        gemmi::Residue* receiver2 = nullptr;
        
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                for (auto& res : chain.residues) {
                    if (res.seqid.num.value == 12 && chain.name == "A") {
                        provider2 = &res;
                    }
                    if (res.seqid.num.value == 13 && chain.name == "A") {
                        receiver2 = &res;
                    }
                }
            }
        }
        
        ASSERT_NE(provider2, nullptr) << "Could not find provider2 residue (A/12)";
        ASSERT_NE(receiver2, nullptr) << "Could not find receiver2 residue (A/13)";
        
        // Clear provider segid to empty
        provider2->segment = "";
        
        bool result = coot::copy_segid_gemmi(*provider2, *receiver2);
        EXPECT_TRUE(result);
        
        // Verify empty segid was copied
        bool empty_segid_copied = receiver2->segment.empty();
        EXPECT_TRUE(empty_segid_copied);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}