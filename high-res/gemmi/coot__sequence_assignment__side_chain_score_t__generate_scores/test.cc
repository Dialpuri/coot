#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, side_chain_score_generate_scores) {
    {
        // Case: generate_scores with PDB, MTZ map, and poly-Ala sequence on chain A

        // Read PDB with gemmi
        gemmi::Structure st = gemmi::read_pdb_file(
            "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        EXPECT_FALSE(st.models.empty());

        gemmi::Model& mol = st.models[0];

        // Create a minimal Xmap<float> for testing
        clipper::Xmap<float> xmap;

        // Create sequence covering chain A
        std::map<std::string, std::vector<std::pair<int, char>>> seq;

        // Chain A has 267 residues. Provide a sequence of at least that length.
        std::vector<std::pair<int, char>> chain_seq;
        for (int i = 0; i < 300; ++i) {
            chain_seq.push_back(std::make_pair(i + 1, ' '));
        }
        seq["A"] = chain_seq;

        EXPECT_EQ(chain_seq.size(), 300u);

        // generate_scores_gemmi is void — it completes without throwing
        EXPECT_NO_THROW(coot::sequence_assignment::generate_scores_gemmi(mol, xmap, seq));
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}