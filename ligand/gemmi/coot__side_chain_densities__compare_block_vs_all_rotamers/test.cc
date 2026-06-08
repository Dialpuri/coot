#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <map>
#include <string>
#include <vector>
#include <utility>
#include <clipper/core/xmap.h>
#include "function.hh"

TEST(CompareBlockVsAllRotamersGemmiTest, basic) {
  // Load a structure
  gemmi::Structure st = gemmi::read_pdb_file(
      "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

  ASSERT_FALSE(st.models.empty());
  gemmi::Model &model = st.models[0];

  // Find residue ASN at chain A, seq 3
  gemmi::Residue *res_p = nullptr;
  for (gemmi::Chain &chain : model.chains) {
    if (chain.name == "A") {
      for (gemmi::Residue &res : chain.residues) {
        if (res.seqid.num.value == 3 && res.name == "ASN") {
          res_p = &res;
          break;
        }
      }
    }
  }

  ASSERT_NE(res_p, nullptr);

  EXPECT_EQ(res_p->seqid.num.value, 3);
  EXPECT_EQ(res_p->name, "ASN");

  // Create a density_box_t (from the compare_block_vs_rotamer gemmi port)
  coot::density_box_t block;
  block.residue_p = res_p;

  // rotamer_limits.first=false means no limits (do everything)
  std::pair<bool, std::vector<std::pair<std::string, std::string>>> rotamer_limits;
  rotamer_limits.first = false;

  // Use a non-existent data_dir so glob_files returns empty
  // -> the function iterates zero dirs, returns empty map
  std::string data_dir = "/nonexistent/path/for/test";

  // Default-construct an empty Xmap
  clipper::Xmap<float> xmap;

  std::map<std::string, std::pair<std::string, double>> result =
      coot::side_chain_densities::compare_block_vs_all_rotamers_gemmi(
          block, res_p, data_dir, rotamer_limits, xmap);

  // With a non-existent directory, glob_files returns empty, so result is empty
  EXPECT_TRUE(result.empty());
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}