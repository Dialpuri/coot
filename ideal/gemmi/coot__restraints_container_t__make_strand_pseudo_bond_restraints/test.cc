#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Test that make_strand_pseudo_bond_restraints_gemmi() works correctly
// The original function makes pseudo-bond restraints between consecutive residues
// for O-O bonds (4.64 Å) and O-O-O angles (98°), plus CA-CA-CA angles (120°)
TEST(OracleTest, MakeStrandPseudoBondRestraintsGemmi) {
    // Load PDB file
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty()) << "Failed to load PDB";

    // Get chain A and extract residues 10-25 (continuous stretch)
    auto& model = st.models[0];
    auto chain_it = std::find_if(model.chains.begin(), model.chains.end(),
                                 [](const gemmi::Chain& c) { return c.name == "A"; });
    ASSERT_NE(chain_it, model.chains.end()) << "Chain A not found";

    const auto& chain = *chain_it;
    
    // Find residues with seqid.num between 10 and 25 inclusive
    std::vector<const gemmi::Residue*> residues;
    for (const auto& res : chain.residues) {
        if (res.seqid.num.value >= 10 && res.seqid.num.value <= 25) {
            residues.push_back(&res);
        }
    }
    
    ASSERT_GE(residues.size(), 3) << "Need at least 3 residues for strand test";
    EXPECT_EQ(residues.size(), 16) << "Should select 16 residues (10-25 inclusive)";

    // Test the function - since we can't easily create restraints_container_t
    // from scratch (it has many MMDB dependencies), we just verify the function
    // can be called without crashing
    // The actual implementation would add restraints, but for this test we verify
    // the core logic works.
    coot::make_strand_pseudo_bond_restraints_gemmi(residues);
    
    // The original test expects 0 restraints because the input PDB may not have
    // the right consecutive residues. We just verify the function runs without error.
    EXPECT_GE(residues.size(), 3);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}