#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// We need to test the same behavior as the original MMDB test, but using gemmi.
// Since gemmi has no Manager, we rewrite the tests to work with Structure.

TEST(OracleTest, SingleBackupGetMolGemmi) {
   coot::molecule_t::modification_info_t::save_info_t si(
       "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb",
       "initial-state");
   const gemmi::Structure* st = si.get_mol_gemmi();
   EXPECT_NE(st, nullptr);
   if (st) {
      EXPECT_EQ(st->models.size(), 1);
      EXPECT_EQ(st->models[0].chains.size(), 2);
   }
}

TEST(OracleTest, MultipleBackupsGetMolGemmi) {
   std::vector<std::string> labels = {"state-0", "state-1", "state-2"};
   for (size_t i = 0; i < labels.size(); i++) {
      coot::molecule_t::modification_info_t::save_info_t si(
          "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb",
          labels[i]);
      const gemmi::Structure* st = si.get_mol_gemmi();
      EXPECT_NE(st, nullptr);
      if (st) {
         EXPECT_EQ(st->models.size(), 1);
         EXPECT_EQ(st->models[0].chains.size(), 2);
      }
   }
}

TEST(OracleTest, InvalidBackupFileHandlingGemmi) {
   coot::molecule_t::modification_info_t::save_info_t si(
       "/nonexistent/path/to/file.pdb",
       "invalid-state");
   const gemmi::Structure* st = si.get_mol_gemmi();
   EXPECT_EQ(st, nullptr);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}