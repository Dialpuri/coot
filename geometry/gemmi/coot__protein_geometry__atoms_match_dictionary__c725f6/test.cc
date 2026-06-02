#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper to find a residue by chain name and sequence number
static const gemmi::Residue* find_residue(
    const gemmi::Structure& st, const std::string& chain_id, int seq_num) {
  for (const auto& model : st.models) {
    for (const auto& chain : model.chains) {
      if (chain.name == chain_id) {
        for (const auto& res : chain.residues) {
          if (res.seqid.num.value == seq_num) {
            return &res;
          }
        }
      }
    }
  }
  return nullptr;
}

TEST(OracleTestGemmi, atoms_match_dictionary) {
  gemmi::Structure st = gemmi::read_pdb_file(
      "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

  // Case 1: Valid standard amino acid residue (ILE at chain A, residue 10)
  {
    const gemmi::Residue* res = find_residue(st, "A", 10);
    ASSERT_NE(res, nullptr);
    EXPECT_EQ(res->name, "ILE");

    std::pair<bool, std::vector<std::string>> result =
        coot::atoms_match_dictionary_gemmi(*res, false, false);

    EXPECT_EQ(result.first, true);
    EXPECT_EQ(result.second.size(), 8u);
    EXPECT_EQ(result.second[0], " N  ");
    EXPECT_EQ(result.second[1], " CA ");
    EXPECT_EQ(result.second[2], " C  ");
    EXPECT_EQ(result.second[3], " O  ");
    EXPECT_EQ(result.second[4], " CB ");
    EXPECT_EQ(result.second[5], " CG1");
    EXPECT_EQ(result.second[6], " CG2");
    EXPECT_EQ(result.second[7], " CD1");
  }

  // Case 2: Residue that does not exist (edge case)
  {
    const gemmi::Residue* res = find_residue(st, "A", 9999);
    EXPECT_EQ(res, nullptr);
  }

  // Case 3: Check hydrogens too
  {
    const gemmi::Residue* res = find_residue(st, "A", 10);
    ASSERT_NE(res, nullptr);
    EXPECT_EQ(res->name, "ILE");

    std::pair<bool, std::vector<std::string>> result =
        coot::atoms_match_dictionary_gemmi(*res, true, true);

    EXPECT_EQ(result.first, true);
    EXPECT_EQ(result.second.size(), 19u);
    EXPECT_EQ(result.second[0], " N  ");
    EXPECT_EQ(result.second[1], " CA ");
    EXPECT_EQ(result.second[2], " C  ");
    EXPECT_EQ(result.second[3], " O  ");
    EXPECT_EQ(result.second[4], " CB ");
    EXPECT_EQ(result.second[5], " CG1");
    EXPECT_EQ(result.second[6], " CG2");
    EXPECT_EQ(result.second[7], " CD1");
    EXPECT_EQ(result.second[8], " H  ");
    EXPECT_EQ(result.second[9], " HA ");
    EXPECT_EQ(result.second[10], " HB ");
    EXPECT_EQ(result.second[11], "HG12");
    EXPECT_EQ(result.second[12], "HG13");
    EXPECT_EQ(result.second[13], "HG21");
    EXPECT_EQ(result.second[14], "HG22");
    EXPECT_EQ(result.second[15], "HG23");
    EXPECT_EQ(result.second[16], "HD11");
    EXPECT_EQ(result.second[17], "HD12");
    EXPECT_EQ(result.second[18], "HD13");
  }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}