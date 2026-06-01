#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, alt_confs_in_molecule) {
  // Case 1: Load example.pdb (standard protein, no explicit alt locs)
  {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    std::vector<std::string> alts = coot::util::alt_confs_in_molecule_gemmi(st);

    EXPECT_EQ(alts.size(), 1u);
    EXPECT_EQ(alts[0], "");
  }

  // Case 2: Load PDB with explicit hydrogen atoms (may have alt locs)
  {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    std::vector<std::string> alts = coot::util::alt_confs_in_molecule_gemmi(st);

    EXPECT_EQ(alts.size(), 1u);
    EXPECT_EQ(alts[0], "");
  }

  // Case 3: Empty/invalid file (mol with no models)
  {
    gemmi::Structure st;
    // Don't load anything — no models

    std::vector<std::string> alts = coot::util::alt_confs_in_molecule_gemmi(st);

    EXPECT_EQ(alts.size(), 0u);
  }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}