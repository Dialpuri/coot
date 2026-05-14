#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, insert_coords_internal_with_new_chain) {
   gemmi::Structure target_st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
   
   const gemmi::Model& target_model_const = target_st.models[0];
   int initial_chains = target_model_const.chains.size();
   int initial_residues = 0;
   for (const auto& chain : target_model_const.chains) {
      initial_residues += chain.residues.size();
   }
   
   gemmi::Structure src_st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
   gemmi::Model& src_model = src_st.models[0];
   if (!src_model.chains.empty()) {
      src_model.chains[0].name = "X";
   }
   
   bool result = coot::insert_coords_internal_gemmi(src_model, target_st.models[0]);
   
   int final_chains = target_st.models[0].chains.size();
   int final_residues = 0;
   int new_chain_residues = 0;
   for (const auto& chain : target_st.models[0].chains) {
      final_residues += chain.residues.size();
      if (chain.name == "X") {
         new_chain_residues = chain.residues.size();
      }
   }
   
   EXPECT_EQ(initial_chains, 2);
   EXPECT_EQ(initial_residues, 268);
   EXPECT_EQ(final_chains, 3);
   EXPECT_EQ(final_residues, 536);
   EXPECT_EQ(new_chain_residues, 267);
   EXPECT_TRUE(final_chains > initial_chains);
   EXPECT_TRUE(final_residues > initial_residues);
}

int main(int argc, char** argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}