#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(PairResidueAtomsTest, adjacent_residues_match_backbone) {
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_GE(st.models.size(), 1) << "no models loaded";

    // Chain 0, residue 0 and 1
    auto &chain = st.models[0].chains[0];
    ASSERT_GE(chain.residues.size(), 2) << "need at least 2 residues";

    const gemmi::Residue &r1 = chain.residues[0];
    const gemmi::Residue &r2 = chain.residues[1];

    std::vector<std::pair<int, int>> pairs =
        coot::util::pair_residue_atoms_gemmi(&r1, &r2);

    EXPECT_EQ(pairs.size(), 6);
    EXPECT_EQ(pairs[0], std::make_pair(0, 0));
    EXPECT_EQ(pairs[1], std::make_pair(1, 1));
    EXPECT_EQ(pairs[2], std::make_pair(2, 2));
    EXPECT_EQ(pairs[3], std::make_pair(3, 3));
    EXPECT_EQ(pairs[4], std::make_pair(4, 4));
    EXPECT_EQ(pairs[5], std::make_pair(5, 5));
}

TEST(PairResidueAtomsTest, same_residue_self_match) {
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_GE(st.models.size(), 1) << "no models loaded";

    auto &chain = st.models[0].chains[0];
    ASSERT_GE(chain.residues.size(), 1) << "need at least 1 residue";

    const gemmi::Residue &r = chain.residues[0];

    std::vector<std::pair<int, int>> pairs =
        coot::util::pair_residue_atoms_gemmi(&r, &r);

    EXPECT_EQ(pairs.size(), 8);
    EXPECT_EQ(pairs[0], std::make_pair(0, 0));
    EXPECT_EQ(pairs[1], std::make_pair(1, 1));
    EXPECT_EQ(pairs[2], std::make_pair(2, 2));
    EXPECT_EQ(pairs[3], std::make_pair(3, 3));
    EXPECT_EQ(pairs[4], std::make_pair(4, 4));
    EXPECT_EQ(pairs[5], std::make_pair(5, 5));
    EXPECT_EQ(pairs[6], std::make_pair(6, 6));
    EXPECT_EQ(pairs[7], std::make_pair(7, 7));
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}