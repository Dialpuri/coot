#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, progressive_residues_in_chain_check) {
  gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

  // gemmi models are 0-indexed; original used GetModel(1) which is model index 0
  gemmi::Model* model = &st.models[0];
  ASSERT_NE(model, nullptr);

  // Case 1: Chain A
  {
    const gemmi::Chain* chain_a = nullptr;
    for (const auto& ch : model->chains) {
      if (ch.name == "A") {
        chain_a = &ch;
        break;
      }
    }
    ASSERT_NE(chain_a, nullptr);
    EXPECT_EQ(chain_a->name, std::string("A"));
    EXPECT_EQ(chain_a->residues.size(), 267);
    short int result = coot::progressive_residues_in_chain_check_gemmi(*chain_a);
    EXPECT_TRUE(result);
  }

  // Case 2: Chain B
  {
    const gemmi::Chain* chain_b = nullptr;
    for (const auto& ch : model->chains) {
      if (ch.name == "B") {
        chain_b = &ch;
        break;
      }
    }
    ASSERT_NE(chain_b, nullptr);
    EXPECT_EQ(chain_b->name, std::string("B"));
    EXPECT_EQ(chain_b->residues.size(), 1);
    short int result = coot::progressive_residues_in_chain_check_gemmi(*chain_b);
    EXPECT_TRUE(result);
  }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}