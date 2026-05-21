#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, CompareResidues) {
    // Read the test PDB file
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Get three residues for testing by finding them in the first model
    const gemmi::Model& model = st.models[0];
    
    const gemmi::Residue* res1 = nullptr;
    const gemmi::Residue* res2 = nullptr;
    const gemmi::Residue* res3 = nullptr;
    
    for (const gemmi::Chain& chain : model.chains) {
        for (const gemmi::Residue& res : chain.residues) {
            if (res.seqid.num.value == 10) res1 = &res;
            if (res.seqid.num.value == 20) res2 = &res;
            if (res.seqid.num.value == 30) res3 = &res;
        }
    }
    
    ASSERT_NE(res1, nullptr) << "Failed to get residue seqnum 10";
    ASSERT_NE(res2, nullptr) << "Failed to get residue seqnum 20";
    ASSERT_NE(res3, nullptr) << "Failed to get residue seqnum 30";
    
    // Test case 1: res1 (seqnum 10) < res2 (seqnum 20) → should return true
    {
        std::pair<const gemmi::Residue*, int> a(res1, 0);
        std::pair<const gemmi::Residue*, int> b(res2, 0);
        bool result = coot::util::compare_residues_gemmi(a, b);
        EXPECT_TRUE(result) << "res1 (seqnum 10) should be less than res2 (seqnum 20)";
    }
    
    // Test case 2: res2 (seqnum 20) < res3 (seqnum 30) → should return true
    {
        std::pair<const gemmi::Residue*, int> a(res2, 0);
        std::pair<const gemmi::Residue*, int> b(res3, 0);
        bool result = coot::util::compare_residues_gemmi(a, b);
        EXPECT_TRUE(result) << "res2 (seqnum 20) should be less than res3 (seqnum 30)";
    }
    
    // Test case 3: res3 (seqnum 30) > res1 (seqnum 10) → should return false
    {
        std::pair<const gemmi::Residue*, int> a(res3, 0);
        std::pair<const gemmi::Residue*, int> b(res1, 0);
        bool result = coot::util::compare_residues_gemmi(a, b);
        EXPECT_FALSE(result) << "res3 (seqnum 30) should not be less than res1 (seqnum 10)";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}