#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__molecule_t__cid_to_residue/gemmi/function.hh"

// Helper: find chain name for a given residue pointer in a model
static std::string get_chain_name(const gemmi::Model& model, gemmi::Residue* res) {
  for (const auto& chain : model.chains) {
    for (const auto& r : chain.residues) {
      if (&r == res) return chain.name;
    }
  }
  return "";
}

TEST(OracleTest, residue_run_t_add_residue) {
  gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
  gemmi::Model& model = st.models[0];

  // Get residue via cid_to_residue_gemmi
  gemmi::Residue* res1 = coot::cid_to_residue_gemmi("//A/1", model);
  ASSERT_NE(res1, nullptr);

  // Verify inputs match oracle
  EXPECT_EQ(get_chain_name(model, res1), std::string("A"));
  EXPECT_EQ(res1->seqid.num.value, 1);
  EXPECT_EQ(res1->name, std::string("MET"));

  // Case 1: add one valid residue
  {
    coot::util::map_to_model_correlation_stats_per_residue_run::residue_run_t run;
    EXPECT_EQ(run.residues.size(), 0u);

    run.add_residue_gemmi(res1);

    EXPECT_EQ(run.residues.size(), 1u);
  }

  // Case 2: add two residues to same run
  {
    coot::util::map_to_model_correlation_stats_per_residue_run::residue_run_t run;
    gemmi::Residue* res2 = coot::cid_to_residue_gemmi("//B/1", model);

    EXPECT_EQ(res1->seqid.num.value, 1);
    EXPECT_EQ(res2 ? res2->seqid.num.value : -1, 1);
    EXPECT_EQ(run.residues.size(), 0u);

    run.add_residue_gemmi(res1);
    EXPECT_EQ(run.residues.size(), 1u);

    run.add_residue_gemmi(res2);
    EXPECT_EQ(run.residues.size(), 2u);
  }

  // Case 3: add nullptr (edge case)
  {
    coot::util::map_to_model_correlation_stats_per_residue_run::residue_run_t run;
    EXPECT_EQ(run.residues.size(), 0u);

    run.add_residue_gemmi(nullptr);

    EXPECT_EQ(run.residues.size(), 1u);
  }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}