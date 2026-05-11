#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper function to find a residue by chain name and residue number
static gemmi::CRA find_residue(gemmi::Model& model, const std::string& chain_name, int seq_num) {
    for (auto& chain : model.chains) {
        if (chain.name == chain_name) {
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == seq_num) {
                    return {&chain, &res, nullptr};
                }
            }
        }
    }
    return {nullptr, nullptr, nullptr};
}

TEST(OracleTest, ResidueAccessAndBonding) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty()) << "PDB loading failed";

    // Create an atom_overlaps_container_t
    coot::atom_overlaps_container_t aoc(&st.models[0], nullptr, false, 0.5, 1.0);

    // Find residues using gemmi API
    gemmi::CRA cra_1 = find_residue(st.models[0], "A", 10);
    gemmi::CRA cra_2 = find_residue(st.models[0], "A", 11);
    gemmi::CRA cra_3 = find_residue(st.models[0], "A", 99);

    // case: cra_1 //A/10
    ASSERT_NE(cra_1.residue, nullptr);
    ASSERT_NE(cra_1.chain, nullptr);
    EXPECT_STREQ(cra_1.chain->name.c_str(), "A");
    EXPECT_EQ(cra_1.residue->seqid.num.value, 10);
    EXPECT_STREQ(cra_1.residue->name.c_str(), "ILE");

    // case: cra_2 //A/11
    ASSERT_NE(cra_2.residue, nullptr);
    ASSERT_NE(cra_2.chain, nullptr);
    EXPECT_STREQ(cra_2.chain->name.c_str(), "A");
    EXPECT_EQ(cra_2.residue->seqid.num.value, 11);
    EXPECT_STREQ(cra_2.residue->name.c_str(), "GLY");

    // case: cra_3 //A/99
    ASSERT_NE(cra_3.residue, nullptr);
    ASSERT_NE(cra_3.chain, nullptr);
    EXPECT_STREQ(cra_3.chain->name.c_str(), "A");
    EXPECT_EQ(cra_3.residue->seqid.num.value, 99);
    EXPECT_STREQ(cra_3.residue->name.c_str(), "ILE");

    // Call are_bonded_residues_gemmi
    bool bonded_1_2 = aoc.are_bonded_residues_gemmi(cra_1, cra_2);
    bool bonded_1_3 = aoc.are_bonded_residues_gemmi(cra_1, cra_3);

    EXPECT_TRUE(bonded_1_2);
    EXPECT_FALSE(bonded_1_3);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}