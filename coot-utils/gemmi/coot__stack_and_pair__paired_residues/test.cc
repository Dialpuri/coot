#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include <utility>
#include <vector>

TEST(OracleTest, stack_and_pair_paired_residues) {
    // Setup: load PDB and initialize geometry
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
        ASSERT_GE(st.models.size(), 1u);

        coot::protein_geometry geom;
        geom.init_standard();

        int imol_enc = 0;

        // Get model 0 (gemmi is 0-indexed, original was model 1 -> index 0)
        gemmi::Model &model = st.models[0];

        // Build residues_vec from all chains/residues in model 0
        std::vector<std::pair<bool, gemmi::Residue *>> residues_vec;
        for (auto &chain : model.chains) {
            for (auto &res : chain.residues) {
                // bool = is_moving (true = moving)
                residues_vec.push_back(std::make_pair(true, &res));
            }
        }

        // Case 1: Call with protein residues (no nucleic acids expected)
        {
            EXPECT_EQ(residues_vec.size(), 268u);
            EXPECT_EQ(imol_enc, 0);

            std::vector<coot::stack_and_pair::paired_residues_info_t> result =
                coot::stack_and_pair::paired_residues_gemmi(model, residues_vec, true, geom, imol_enc);

            EXPECT_EQ(result.size(), 0u);
        }

        // Case 2: Call with empty residues_vec
        {
            std::vector<std::pair<bool, gemmi::Residue *>> empty_vec;
            EXPECT_EQ(empty_vec.size(), 0u);

            std::vector<coot::stack_and_pair::paired_residues_info_t> result =
                coot::stack_and_pair::paired_residues_gemmi(model, empty_vec, false, geom, imol_enc);

            EXPECT_EQ(result.size(), 0u);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}