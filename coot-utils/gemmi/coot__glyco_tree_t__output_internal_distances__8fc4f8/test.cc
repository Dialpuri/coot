#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "function.hh"

// Find a residue by chain ID and sequence number
static const gemmi::Residue* find_residue(const gemmi::Model &model, const std::string &chain_id, int seq_num) {
  for (const auto &chain : model.chains) {
    if (chain.name == chain_id) {
      for (const auto &res : chain.residues) {
        if (res.seqid.num.value == seq_num) {
          return &res;
        }
      }
    }
  }
  return nullptr;
}

TEST(GemmiTest, glyco_tree_output_internal_distances) {
  gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
  ASSERT_FALSE(st.models.empty());

  gemmi::Model &model = st.models[0];

  // Get two adjacent residues
  const gemmi::Residue *res1_p = find_residue(model, "A", 1);
  const gemmi::Residue *res2_p = find_residue(model, "A", 2);

  ASSERT_NE(res1_p, nullptr);
  ASSERT_NE(res2_p, nullptr);

  double dist_crit = 8.0;
  int model_number = 1;
  std::string chain_id = "A";

  std::vector<const gemmi::Residue *> residues_list;
  residues_list.push_back(res1_p);
  residues_list.push_back(res2_p);

  // Case 1: two adjacent residues, moderate cutoff
  {
    std::string tmpfile = "/tmp/glyco_dist_test_1.txt";
    {
      std::ofstream f(tmpfile.c_str());
      ASSERT_TRUE(f.is_open());
      coot::output_internal_distances_gemmi(*res1_p, residues_list, dist_crit, f, model_number, chain_id);
    }

    // Read back the file contents and check first 5 lines and total count
    std::ifstream rf(tmpfile.c_str());
    std::string line;
    int line_count = 0;
    std::vector<std::string> first_lines;
    while (std::getline(rf, line)) {
      line_count++;
      if (line_count <= 5) {
        first_lines.push_back(line);
      }
    }

    EXPECT_EQ(line_count, 59);

    EXPECT_FALSE(first_lines[0].empty());
    EXPECT_NE(first_lines[0].find("[spec: model 1 \"A\"    2 \"\" \" N  \" \"\"]"), std::string::npos);
    EXPECT_NE(first_lines[0].find("2.73041"), std::string::npos);

    EXPECT_NE(first_lines[1].find("[spec: model 1 \"A\"    2 \"\" \" N  \" \"\"]"), std::string::npos);
    EXPECT_NE(first_lines[1].find("3.24333"), std::string::npos);

    EXPECT_NE(first_lines[2].find("3.8512"), std::string::npos);

    EXPECT_NE(first_lines[3].find("3.76242"), std::string::npos);

    EXPECT_NE(first_lines[4].find("4.2848"), std::string::npos);
  }

  // Case 2: very tight cutoff — should produce few or no output
  {
    std::string tmpfile = "/tmp/glyco_dist_test_2.txt";
    double tight = 2.70;
    {
      std::ofstream f(tmpfile.c_str());
      ASSERT_TRUE(f.is_open());
      coot::output_internal_distances_gemmi(*res1_p, residues_list, tight, f, model_number, chain_id);
    }

    std::ifstream rf(tmpfile.c_str());
    std::string line;
    int line_count = 0;
    while (std::getline(rf, line)) {
      line_count++;
    }
    EXPECT_EQ(line_count, 0);
  }

  // Case 3: single residue in list (self-skip) — should produce 0 lines
  {
    std::string tmpfile = "/tmp/glyco_dist_test_3.txt";
    std::vector<const gemmi::Residue *> single_list;
    single_list.push_back(res1_p);
    {
      std::ofstream f(tmpfile.c_str());
      ASSERT_TRUE(f.is_open());
      coot::output_internal_distances_gemmi(*res1_p, single_list, dist_crit, f, model_number, chain_id);
    }

    std::ifstream rf(tmpfile.c_str());
    std::string line;
    int line_count = 0;
    while (std::getline(rf, line)) {
      line_count++;
    }
    EXPECT_EQ(line_count, 0);
  }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}