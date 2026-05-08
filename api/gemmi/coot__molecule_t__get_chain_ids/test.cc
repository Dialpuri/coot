#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, GetChainIds) {
   // Test 1: example.pdb (no hydrogens)
   {
      gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
      EXPECT_EQ(st.models.size(), 1);
      EXPECT_EQ(st.models[0].chains.size(), 2);
      EXPECT_EQ(st.models[0].chains[0].name, "A");
      EXPECT_EQ(st.models[0].chains[1].name, "B");
   }
   
   // Test 2: example-hydrogen.pdb (with hydrogens)
   {
      gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
      EXPECT_EQ(st.models.size(), 1);
      EXPECT_EQ(st.models[0].chains.size(), 2);
      EXPECT_EQ(st.models[0].chains[0].name, "A");
      EXPECT_EQ(st.models[0].chains[1].name, "B");
   }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
