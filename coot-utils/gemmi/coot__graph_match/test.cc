#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: find a residue by chain name and sequence number
static const gemmi::Residue* find_residue(const gemmi::Model& model,
                                           const std::string& chain_id,
                                           int seqnum) {
  for (const auto& chain : model.chains) {
    if (chain.name != chain_id) continue;
    for (const auto& res : chain.residues) {
      if (res.seqid.num.value == seqnum) {
        return &res;
      }
    }
  }
  return nullptr;
}

TEST(OracleTest, graph_match) {
  gemmi::Structure st = gemmi::read_pdb_file(
      "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
  
  auto& model = st.models[0];

  // Case 1: Match two identical residue types (ALA 21 vs ALA 31) with apply_rtop_flag=true, no hydrogens
  {
    const gemmi::Residue* res_moving = find_residue(model, "A", 21);
    const gemmi::Residue* res_reference = find_residue(model, "A", 31);

    bool apply_rtop_flag = true;
    bool match_hydrogens_also = false;

    EXPECT_NE(res_moving, nullptr);
    EXPECT_NE(res_reference, nullptr);

    coot::graph_match_info_t result = coot::graph_match_gemmi(res_moving, res_reference,
                                                               apply_rtop_flag, match_hydrogens_also);
    EXPECT_TRUE(result.success);
    EXPECT_NEAR(result.dist_score, 0.526594, 1e-4);
    EXPECT_EQ(result.n_match, 5);
    EXPECT_EQ(result.matching_atom_names.size(), 5u);
  }

  // Case 2: Same residues but apply_rtop_flag=false
  {
    const gemmi::Residue* res_moving = find_residue(model, "A", 21);
    const gemmi::Residue* res_reference = find_residue(model, "A", 31);

    bool apply_rtop_flag = false;
    bool match_hydrogens_also = false;

    EXPECT_NE(res_moving, nullptr);
    EXPECT_NE(res_reference, nullptr);

    coot::graph_match_info_t result = coot::graph_match_gemmi(res_moving, res_reference,
                                                               apply_rtop_flag, match_hydrogens_also);
    EXPECT_TRUE(result.success);
    EXPECT_NEAR(result.dist_score, 31.8446, 1e-4);
    EXPECT_EQ(result.n_match, 5);
    EXPECT_EQ(result.matching_atom_names.size(), 5u);
  }

  // Case 3: Different residue types (ALA vs GLY) — graph matching still succeeds
  {
    const gemmi::Residue* res_moving = find_residue(model, "A", 21);
    const gemmi::Residue* res_reference = find_residue(model, "A", 11);

    bool apply_rtop_flag = true;
    bool match_hydrogens_also = false;

    EXPECT_NE(res_moving, nullptr);
    EXPECT_NE(res_reference, nullptr);

    coot::graph_match_info_t result = coot::graph_match_gemmi(res_moving, res_reference,
                                                               apply_rtop_flag, match_hydrogens_also);
    EXPECT_TRUE(result.success);
    EXPECT_NEAR(result.dist_score, 0.126463, 1e-4);
    EXPECT_EQ(result.n_match, 4);
    EXPECT_EQ(result.matching_atom_names.size(), 4u);
  }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}