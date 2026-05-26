#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, GetPhiPsi) {
    // Read example PDB file
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    ASSERT_FALSE(st.models.empty()) << "Could not get model";
    
    gemmi::Model& model = st.models[0];
    
    // Get chain A
    bool found_chain = false;
    gemmi::Chain* chain = nullptr;
    for (gemmi::Chain& c : model.chains) {
        if (c.name == "A") {
            chain = &c;
            found_chain = true;
            break;
        }
    }
    ASSERT_TRUE(found_chain) << "Could not get chain A";
    
    // Access residues (gemmi uses 0-indexed, but example.pdb residues should be consecutive)
    ASSERT_GE(chain->residues.size(), 3) << "Need at least 3 residues";
    
    const gemmi::Residue* residues[3];
    residues[0] = &chain->residues[1];  // GLU (seqid.num.value = 2)
    residues[1] = &chain->residues[2];  // ASN (seqid.num.value = 3)
    residues[2] = &chain->residues[3];  // PHE (seqid.num.value = 4)
    
    ASSERT_NE(residues[0], nullptr);
    ASSERT_NE(residues[1], nullptr);
    ASSERT_NE(residues[2], nullptr);
    
    // Test the function that takes 3 residues
    std::pair<bool, coot::util::phi_psi_with_residues_t> result = 
        coot::util::get_phi_psi_gemmi(residues[0], residues[1], residues[2]);
    
    EXPECT_TRUE(result.first) << "get_phi_psi should succeed";
    
    if (result.first) {
        EXPECT_NEAR(result.second.phi_, -75.3333, 1e-4);
        EXPECT_NEAR(result.second.psi_, 1.75129, 1e-4);
        
        // Check residue identifiers
        EXPECT_EQ(residues[0]->seqid.num.value, 2);
        EXPECT_EQ(residues[0]->name, "GLU");
        
        EXPECT_EQ(residues[1]->seqid.num.value, 3);
        EXPECT_EQ(residues[1]->name, "ASN");
        
        EXPECT_EQ(residues[2]->seqid.num.value, 4);
        EXPECT_EQ(residues[2]->name, "PHE");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}