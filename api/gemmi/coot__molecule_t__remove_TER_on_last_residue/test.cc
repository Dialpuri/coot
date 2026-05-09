#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper function to count TER atoms in a residue
int count_TER_atoms(const gemmi::Residue& residue) {
    int count = 0;
    for (const auto& atom : residue.atoms) {
        if (atom.name.find("TER") != std::string::npos) {
            count++;
        }
    }
    return count;
}

TEST(OracleTest, RemoveTEROnLastResidue) {
   // Load a PDB file
   gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
   
   ASSERT_GE(st.models.size(), 1) << "Failed to load PDB file or no models found";
   
   // Get the first model
   gemmi::Model& model = st.models[0];
   
   int nChains = model.chains.size();
   EXPECT_EQ(nChains, 2);
   
   // Test both chains (A and B)
   for (size_t ichain = 0; ichain < model.chains.size(); ichain++) {
      gemmi::Chain& chain = model.chains[ichain];
      
      if (chain.residues.empty()) continue; // Skip if no residues
      
      // Get the last residue
      gemmi::Residue& last_residue = chain.residues.back();
      
      // Count TER atoms in the last residue before
      int n_TER_before = count_TER_atoms(last_residue);
      
      // Call the function under test
      coot::molecule_t::remove_TER_on_last_residue_gemmi(chain);
      
      // Count TER atoms in the last residue after
      int n_TER_after = count_TER_atoms(last_residue);
      
      int n_TER_removed = n_TER_before - n_TER_after;
      
      // Verify the function doesn't crash and the counts are as expected
      EXPECT_EQ(n_TER_before, 0) << "Chain " << ichain << " should have 0 TER atoms before";
      EXPECT_EQ(n_TER_after, 0) << "Chain " << ichain << " should have 0 TER atoms after";
      EXPECT_EQ(n_TER_removed, 0) << "Chain " << ichain << " should have 0 TER atoms removed";
   }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}