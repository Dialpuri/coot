#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_pucker_analysis_info) {
   std::string pdb_path = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";
   gemmi::Structure st = gemmi::read_pdb_file(pdb_path);
   
   // Verify structure loaded correctly
   ASSERT_FALSE(st.models.empty());
   ASSERT_FALSE(st.models[0].chains.empty());
   
   // Count atoms in each residue
   bool has_qualifying_residue = false;
   for (const auto& model : st.models) {
      for (const auto& chain : model.chains) {
         for (const auto& res : chain.residues) {
            if (res.atoms.size() > 14) {
               has_qualifying_residue = true;
               break;
            }
         }
      }
   }
   
   // Test with the structure
   coot::molecule_t mol(st);
   std::string result = mol.get_pucker_analysis_info_gemmi();
   
   // Expected: empty string when no residues have >14 atoms
   EXPECT_EQ(result.length(), 0);
}

int main(int argc, char** argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}