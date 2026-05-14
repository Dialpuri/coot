#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, MutateResiduesGemmi) {
  gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

  // Test 1: Mutate ILE to ALA at chain A, residue 10
  {
    // Find residue at A/10
    gemmi::Residue *res10 = nullptr;
    for (auto& model : st.models)
      for (auto& chain : model.chains)
        if (chain.name == "A")
          for (auto& res : chain.residues)
            if (res.seqid.num.value == 10 && res.seqid.icode == ' ') {
              res10 = &res;
              break;
            }
    ASSERT_TRUE(res10 != nullptr);
    std::string current_type = res10->name;

    int status = coot::molecule_t::mutate_gemmi(st, "A", 10, ' ', "ALA");
    EXPECT_EQ(current_type, "ILE");
    EXPECT_EQ(res10->name, "ALA");
    EXPECT_EQ(status, 1);
  }

  // Test 2: Mutate ALA to GLY at same residue (A/10)
  {
    gemmi::Residue *res10 = nullptr;
    for (auto& model : st.models)
      for (auto& chain : model.chains)
        if (chain.name == "A")
          for (auto& res : chain.residues)
            if (res.seqid.num.value == 10 && res.seqid.icode == ' ') {
              res10 = &res;
              break;
            }
    ASSERT_TRUE(res10 != nullptr);
    std::string current_type = res10->name;  // should be "ALA"

    int status = coot::molecule_t::mutate_gemmi(st, "A", 10, ' ', "GLY");
    EXPECT_EQ(current_type, "ALA");
    EXPECT_EQ(res10->name, "GLY");
    EXPECT_EQ(status, 1);
  }

  // Test 3: Invalid residue type (should fail)
  {
    gemmi::Residue *res10 = nullptr;
    for (auto& model : st.models)
      for (auto& chain : model.chains)
        if (chain.name == "A")
          for (auto& res : chain.residues)
            if (res.seqid.num.value == 10 && res.seqid.icode == ' ') {
              res10 = &res;
              break;
            }
    ASSERT_TRUE(res10 != nullptr);
    std::string current_type = res10->name;  // should be "GLY"

    int status = coot::molecule_t::mutate_gemmi(st, "A", 10, ' ', "XXX");
    EXPECT_EQ(current_type, "GLY");  // should stay unchanged
    EXPECT_EQ(res10->name, "GLY");
    EXPECT_EQ(status, -1);
  }

  // Test 4: Mutate residue A/2 from GLU to SER
  {
    gemmi::Residue *res2 = nullptr;
    for (auto& model : st.models)
      for (auto& chain : model.chains)
        if (chain.name == "A")
          for (auto& res : chain.residues)
            if (res.seqid.num.value == 2 && res.seqid.icode == ' ') {
              res2 = &res;
              break;
            }
    ASSERT_TRUE(res2 != nullptr);
    std::string current_type = res2->name;

    int status = coot::molecule_t::mutate_gemmi(st, "A", 2, ' ', "SER");
    EXPECT_EQ(current_type, "GLU");
    EXPECT_EQ(res2->name, "SER");
    EXPECT_EQ(status, 1);
  }
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}