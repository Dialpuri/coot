#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <fstream>
#include <sstream>
#include "function.hh"

TEST(OracleTest, copy_and_add_chain_residues_to_chain) {
   // Load the main PDB (has chains A and B, 267 residues in chain A)
   gemmi::Structure target_st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
   
   // Create a temporary PDB with a single water residue
   std::string temp_pdb_path = "/tmp/test_single_hoh.pdb";
   std::ostringstream oss;
   oss << "ATOM      1  O   HOH X   1       0.000   0.000   0.000  1.00  0.00           O\n";
   oss << "TER       2      HOH X   1\n";
   std::ofstream ofs(temp_pdb_path);
   ofs << oss.str();
   ofs.close();
   
   gemmi::Structure source_st = gemmi::read_pdb_file(temp_pdb_path);

   // Get model 1 (gemmi uses 0-indexing, so model 1 -> index 0)
   ASSERT_FALSE(target_st.models.empty());
   ASSERT_FALSE(source_st.models.empty());
   gemmi::Model& target_model = target_st.models[0];
   gemmi::Model& source_model = source_st.models[0];

   // Get chain A from target (267 residues) and chain X from source (1 water)
   auto target_chain_it = std::find_if(target_model.chains.begin(), target_model.chains.end(),
                                       [](const gemmi::Chain& c) { return c.name == "A"; });
   auto source_chain_it = std::find_if(source_model.chains.begin(), source_model.chains.end(),
                                       [](const gemmi::Chain& c) { return c.name == "X"; });

   ASSERT_NE(target_chain_it, target_model.chains.end());
   ASSERT_NE(source_chain_it, source_model.chains.end());

   gemmi::Chain& target_chain = *target_chain_it;
   gemmi::Chain& source_chain = *source_chain_it;

   int nres_before = static_cast<int>(target_chain.residues.size());
   EXPECT_EQ(nres_before, 267);

   // Call the ported function
   coot::molecule_t::copy_and_add_chain_residues_to_chain_gemmi(source_st, target_st, "X", "A");

   int nres_after = static_cast<int>(target_chain.residues.size());
   EXPECT_EQ(nres_after, 268);

   // Verify the new residue is a water
   gemmi::Residue& new_res = target_chain.residues.back();
   EXPECT_EQ(new_res.name, "HOH");
}

int main(int argc, char** argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}