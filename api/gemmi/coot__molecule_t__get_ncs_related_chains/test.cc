#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_ncs_related_chains) {
   gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
   
   std::vector<std::vector<std::string>> result = coot::molecule_t::get_ncs_related_chains_gemmi(st);
   
   EXPECT_EQ(result.size(), 2u);
   
   {
      // Group 0
      ASSERT_EQ(result[0].size(), 1u);
      EXPECT_EQ(result[0][0], "A");
   }
   
   {
      // Group 1
      ASSERT_EQ(result[1].size(), 1u);
      EXPECT_EQ(result[1][0], "B");
   }
}

int main(int argc, char** argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}