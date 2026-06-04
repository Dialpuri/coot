#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

static gemmi::Residue* find_residue(gemmi::Model& model,
                                     char chain_id,
                                     int seqnum) {
  for (auto& chain : model.chains) {
    if (chain.name.size() && chain.name[0] == chain_id) {
      for (auto& res : chain.residues) {
        if (res.seqid.num.value == seqnum) {
          return &res;
        }
      }
    }
  }
  return nullptr;
}

TEST(OracleTestGemmi, add_xH_H_case1_SER_HG) {
  gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-dev/coot/data/tutorial-modern.pdb");
  ASSERT_FALSE(st.models.empty());
  gemmi::Model& model = st.models[0];

  gemmi::Residue* res = find_residue(model, 'A', 3);
  ASSERT_NE(res, nullptr);

  EXPECT_EQ(res->name, "SER");
  EXPECT_EQ(res->seqid.num.value, 3);

  int atoms_before = static_cast<int>(res->atoms.size());
  EXPECT_EQ(atoms_before, 6);

  std::string H_name = "HG";
  std::string at_name_1 = "OG";
  std::string at_name_2 = "CB";
  std::string at_name_3 = "CA";
  double bl = 0.96;
  double angle = 109.5;
  double tor_initial = 0.0;

  std::vector<gemmi::Atom *> result =
    coot::reduce::add_xH_H_gemmi(H_name, at_name_1, at_name_2, at_name_3,
                                 bl, angle, tor_initial, *res);

  EXPECT_EQ(result.size(), 1u);
  EXPECT_NE(result.empty(), true);
  // Atom name in gemmi may be padded; compare trimmed
  std::string found_name = result[0]->name;
  // Trim whitespace
  size_t fs = found_name.find_first_not_of(" \t");
  size_t fl = found_name.find_last_not_of(" \t");
  if (fs != std::string::npos)
    found_name = found_name.substr(fs, fl - fs + 1);
  EXPECT_EQ(found_name, "HG");

  int atoms_after = static_cast<int>(res->atoms.size());
  EXPECT_EQ(atoms_after, 7);
}

TEST(OracleTestGemmi, add_xH_H_case2_invalid_atoms) {
  gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-dev/coot/data/tutorial-modern.pdb");
  ASSERT_FALSE(st.models.empty());
  gemmi::Model& model = st.models[0];

  gemmi::Residue* res = find_residue(model, 'A', 3);
  ASSERT_NE(res, nullptr);

  EXPECT_EQ(res->name, "SER");
  EXPECT_EQ(res->seqid.num.value, 3);

  std::string H_name = "HX";
  std::string at_name_1 = "NONEXISTENT";
  std::string at_name_2 = "ALSO_MISSING";
  std::string at_name_3 = "AND_THIS";
  double bl = 0.96;
  double angle = 109.5;
  double tor_initial = 0.0;

  std::vector<gemmi::Atom *> result =
    coot::reduce::add_xH_H_gemmi(H_name, at_name_1, at_name_2, at_name_3,
                                 bl, angle, tor_initial, *res);

  EXPECT_EQ(result.size(), 0u);
  EXPECT_TRUE(result.empty());
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}