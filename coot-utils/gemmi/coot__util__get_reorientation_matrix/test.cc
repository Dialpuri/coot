#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: find a residue by chain name and sequence number
static const gemmi::Residue* find_residue(const gemmi::Model& model,
                                           const std::string& chain_name,
                                           int seqnum) {
  for (const auto& chain : model.chains) {
    if (chain.name == chain_name) {
      for (const auto& res : chain.residues) {
        if (res.seqid.num.value == seqnum) {
          return &res;
        }
      }
    }
  }
  return nullptr;
}

TEST(OracleTest, get_reorientation_matrix) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Use first model
    ASSERT_FALSE(st.models.empty());
    const gemmi::Model& model = st.models[0];

    // Case 1: Valid adjacent residues (10 ILE and 11 GLY in chain A)
    {
        const gemmi::Residue *res_current = find_residue(model, "A", 10);
        const gemmi::Residue *res_next    = find_residue(model, "A", 11);

        EXPECT_NE(res_current, nullptr);
        EXPECT_NE(res_next, nullptr);

        std::pair<bool, clipper::RTop_orth> result =
            coot::util::get_reorientation_matrix_gemmi(res_current, res_next);

        EXPECT_TRUE(result.first);

        std::string rtop_format = result.second.format().data();
        EXPECT_EQ(rtop_format,
"|    0.3238,   -0.4266,    0.8445|\n"
"|    0.7277,    0.6828,    0.0659|\n"
"|   -0.6047,    0.5932,    0.5315|\n"
"(    -35.26,    -18.73,     37.89)");
    }

    // Case 2: Null residue_next (edge case)
    {
        const gemmi::Residue *res_current = find_residue(model, "A", 10);
        const gemmi::Residue *res_next    = nullptr;

        std::pair<bool, clipper::RTop_orth> result =
            coot::util::get_reorientation_matrix_gemmi(res_current, res_next);

        EXPECT_FALSE(result.first);
    }

    // Case 3: Both null (edge case)
    {
        const gemmi::Residue *res_current = nullptr;
        const gemmi::Residue *res_next    = nullptr;

        std::pair<bool, clipper::RTop_orth> result =
            coot::util::get_reorientation_matrix_gemmi(res_current, res_next);

        EXPECT_FALSE(result.first);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}