#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: find a residue by chain name and sequence number (ignoring name/icode)
static gemmi::Residue* find_residue_by_seq(gemmi::Model& model,
                                           const std::string& chain_name,
                                           int seq_num) {
  for (gemmi::Chain& chain : model.chains) {
    if (chain.name != chain_name)
      continue;
    for (gemmi::Residue& res : chain.residues) {
      if (res.seqid.num.value == seq_num)
        return &res;
    }
  }
  return nullptr;
}

TEST(OracleTestGemmi, delete_alt_confs_except) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];

    // Case 1: alt_conf="" (empty) → atoms kept (before == after)
    {
        gemmi::Residue* res = find_residue_by_seq(model, "A", 10);
        ASSERT_NE(res, nullptr);

        int n_before = static_cast<int>(res->atoms.size());
        EXPECT_EQ(n_before, 8);

        std::string alt_conf = "";
        coot::util::delete_alt_confs_except_gemmi(*res, alt_conf);

        int n_after = static_cast<int>(res->atoms.size());
        EXPECT_EQ(n_after, 8);
    }

    // Case 2: alt_conf="A" does NOT match empty altLoc → all atoms deleted
    {
        gemmi::Structure st2 = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        gemmi::Model& model2 = st2.models[0];

        gemmi::Residue* res = find_residue_by_seq(model2, "A", 10);
        ASSERT_NE(res, nullptr);

        int n_before = static_cast<int>(res->atoms.size());
        EXPECT_EQ(n_before, 8);

        std::string alt_conf = "A";
        coot::util::delete_alt_confs_except_gemmi(*res, alt_conf);

        int n_after = static_cast<int>(res->atoms.size());
        EXPECT_EQ(n_after, 0);
    }

    // Case 3: non-existent residue
    {
        gemmi::Residue* res = find_residue_by_seq(model, "X", 9999);
        EXPECT_EQ(res, nullptr);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}