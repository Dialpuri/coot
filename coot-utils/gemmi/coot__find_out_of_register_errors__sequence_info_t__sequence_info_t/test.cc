#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, SequenceInfo) {
   // Load PDB
   gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

   // Get model 0 (gemmi is 0-indexed, original GetModel(1) was 1-indexed)
   ASSERT_FALSE(st.models.empty());
   auto& model = st.models[0];
   ASSERT_FALSE(model.chains.empty());
   auto& chain = model.chains[0];

   // Case 1: construct from residue at index 9 (ILE 10)
   {
      auto& res = chain.residues[9];
      int rn = res.seqid.num.value;
      std::string rn_str = res.name;
      coot::find_out_of_register_errors::sequence_info_t si(&res, rn, rn_str);

      EXPECT_NE(si.residue_p, nullptr);
      EXPECT_EQ(si.res_no, 10);
      EXPECT_EQ(si.res_name, "ILE");
   }

   // Case 2: construct from first residue (MET 1)
   {
      auto& res = chain.residues[0];
      int rn = res.seqid.num.value;
      std::string rn_str = res.name;
      coot::find_out_of_register_errors::sequence_info_t si(&res, rn, rn_str);

      EXPECT_NE(si.residue_p, nullptr);
      EXPECT_EQ(si.res_no, 1);
      EXPECT_EQ(si.res_name, "MET");
   }

   // Case 3: construct with nullptr residue (edge case)
   {
      coot::find_out_of_register_errors::sequence_info_t si(nullptr, 0, "UNK");

      EXPECT_EQ(si.residue_p, nullptr);
      EXPECT_EQ(si.res_no, 0);
      EXPECT_EQ(si.res_name, "UNK");
   }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}