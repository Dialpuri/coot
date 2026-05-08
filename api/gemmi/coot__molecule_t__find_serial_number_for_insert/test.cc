#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, find_serial_number_for_insert) {
   gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
   const gemmi::Model& model = st.models[0];

   // Case 1: insert before residue 10 (should find residue 10)
   {
      std::pair<int, gemmi::CRA> result1 = coot::find_serial_number_for_insert_gemmi(9, "", "A", model);
      EXPECT_EQ(result1.first, 9);
      ASSERT_NE(result1.second.chain, nullptr);
      EXPECT_EQ(result1.second.residue->seqid.num.value, 10);
      EXPECT_EQ(result1.second.residue->seqid.icode, ' ');
   }

   // Case 2: insert between residues (after 10, before 11) with seqnum 10, ins_code 'A'
   {
      std::pair<int, gemmi::CRA> result2 = coot::find_serial_number_for_insert_gemmi(10, "A", "A", model);
      EXPECT_EQ(result2.first, 10);
      ASSERT_NE(result2.second.chain, nullptr);
      EXPECT_EQ(result2.second.residue->seqid.num.value, 11);
      EXPECT_EQ(result2.second.residue->seqid.icode, ' ');
   }

   // Case 3: insert after the last residue (should return -1, empty CRA)
   {
      std::pair<int, gemmi::CRA> result3 = coot::find_serial_number_for_insert_gemmi(999, "", "A", model);
      EXPECT_EQ(result3.first, -1);
      EXPECT_EQ(result3.second.chain, nullptr);
   }

   // Case 4: insert before residue 100 (assuming residue 100 exists)
   {
      std::pair<int, gemmi::CRA> result4 = coot::find_serial_number_for_insert_gemmi(99, "", "A", model);
      EXPECT_EQ(result4.first, 87);
      ASSERT_NE(result4.second.chain, nullptr);
      EXPECT_EQ(result4.second.residue->seqid.num.value, 100);
      EXPECT_EQ(result4.second.residue->seqid.icode, ' ');
   }

   // Case 5: invalid chain
   {
      std::pair<int, gemmi::CRA> result5 = coot::find_serial_number_for_insert_gemmi(10, "", "Z", model);
      EXPECT_EQ(result5.first, -1);
      EXPECT_EQ(result5.second.chain, nullptr);
   }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}