#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "ideal/zo-rama.hh"

TEST(OracleTest, nmer_crankshaft_set) {
  // Setup: load molecule using gemmi
  gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
  ASSERT_FALSE(st.models.empty());
  const gemmi::Model *mol = &st.models[0];

  zo::rama_table_set zorts;
  zorts.init();

  // Case 1: valid residue spec — middle of chain A, residue 10 (SER), n_peptides=3
  {
    coot::residue_spec_t spec("A", 10, "");
    unsigned int n_peptides = 3;

    EXPECT_NO_THROW({
      coot::nmer_crankshaft_set_gemmi cs(spec, n_peptides, zorts, mol);
      EXPECT_EQ(cs.size(), 0u);
    });
  }

  // Case 2: edge — n_peptides=1 (minimum)
  {
    coot::residue_spec_t spec("A", 10, "");
    unsigned int n_peptides = 1;

    EXPECT_NO_THROW({
      coot::nmer_crankshaft_set_gemmi cs(spec, n_peptides, zorts, mol);
      EXPECT_EQ(cs.size(), 0u);
    });
  }

  // Case 3: non-existent chain — should produce empty or exception
  {
    coot::residue_spec_t spec("Z", 10, "");
    unsigned int n_peptides = 3;

    EXPECT_NO_THROW({
      coot::nmer_crankshaft_set_gemmi cs(spec, n_peptides, zorts, mol);
      EXPECT_EQ(cs.size(), 0u);
    });
  }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}