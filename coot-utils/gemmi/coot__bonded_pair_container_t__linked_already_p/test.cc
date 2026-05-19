#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, BondedPairsContainer) {
   // Load PDB to get residues
   gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
   
   // Get first model
   ASSERT_FALSE(st.models.empty()) << "No model found";
   gemmi::Model& model = st.models[0];
   
   // Get chain A (PDB chain ID 'A')
   gemmi::Chain* chain = nullptr;
   for (gemmi::Chain& c : model.chains) {
       if (c.name == "A") {
           chain = &c;
           break;
       }
   }
   ASSERT_NE(chain, nullptr) << "No chain A found";
   
   // Get two residues - find by seqid.num.value
   gemmi::Residue* res1 = nullptr;
   gemmi::Residue* res2 = nullptr;
   for (gemmi::Residue& r : chain->residues) {
       if (r.seqid.num.value == 10) res1 = &r;
       if (r.seqid.num.value == 11) res2 = &r;
   }
   
   ASSERT_NE(res1, nullptr) << "Failed to get residue 10";
   ASSERT_NE(res2, nullptr) << "Failed to get residue 11";
   
   // Create bonded_pair_container_t_gemmi and add a bonded pair
   coot::bonded_pair_container_t_gemmi container;
   
   // First, add a bonded pair between res1 and res2
   {
      coot::bonded_pair_t_gemmi bp(res1, res2, false, false, "test_link");
      // Manually add to vector since try_add doesn't exist in ported version
      container.bonded_residues.push_back(bp);
      
      EXPECT_EQ(container.size(), 1u) << "After adding first pair, container size should be 1";
   }
   
   // Test case 1: check if res1 and res2 are already linked (should be true)
   {
      bool result = container.linked_already_p(res1, res2);
      int res1_seq = res1->seqid.num.value;
      int res2_seq = res2->seqid.num.value;
      EXPECT_TRUE(result) << "res1: \"//A/" << res1_seq << "\" res2: \"//A/" << res2_seq << "\" (already added)";
   }
   
   // Test case 2: check if res1 and res2 in swapped order are linked (should still be true)
   {
      bool result = container.linked_already_p(res2, res1);
      int res1_seq = res1->seqid.num.value;
      int res2_seq = res2->seqid.num.value;
      EXPECT_TRUE(result) << "res1: \"//A/" << res2_seq << "\" res2: \"//A/" << res1_seq << "\" (swapped order)";
   }
   
   // Test case 3: check with a residue not in the container (should be false)
   {
      gemmi::Residue* res3 = nullptr;
      for (gemmi::Residue& r : chain->residues) {
          if (r.seqid.num.value == 20) {
              res3 = &r;
              break;
          }
      }
      int res1_seq = res1->seqid.num.value;
      if (res3) {
         int res3_seq = res3->seqid.num.value;
         bool result = container.linked_already_p(res1, res3);
         EXPECT_FALSE(result) << "res1: \"//A/" << res1_seq << "\" res2: \"//A/" << res3_seq << "\" (not added)";
      }
   }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}