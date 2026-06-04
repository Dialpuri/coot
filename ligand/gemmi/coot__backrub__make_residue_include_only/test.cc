#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, make_residue_include_only) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];

    // Get three consecutive residues: A/9, A/10, A/11
    const gemmi::Chain* chain_a = nullptr;
    for (const auto& ch : model.chains) {
        if (ch.name == "A") { chain_a = &ch; break; }
    }
    ASSERT_NE(chain_a, nullptr);

    const gemmi::Residue* prev_res  = nullptr;  // A/9  LYS
    const gemmi::Residue* this_res  = nullptr;  // A/10 ILE
    const gemmi::Residue* next_res  = nullptr;  // A/11 GLY

    for (const auto& res : chain_a->residues) {
        if (res.seqid.num.value == 9)  prev_res  = &res;
        if (res.seqid.num.value == 10) this_res  = &res;
        if (res.seqid.num.value == 11) next_res  = &res;
    }

    EXPECT_NE(prev_res, nullptr);
    EXPECT_NE(this_res, nullptr);
    EXPECT_NE(next_res, nullptr);

    // Case 1: Valid residue with backbone atoms C and O
    {
        std::vector<std::string> atoms_to_include;
        atoms_to_include.push_back(" C  ");
        atoms_to_include.push_back(" O  ");

        coot::minimol::residue r = coot::backrub::make_residue_include_only_gemmi(prev_res, atoms_to_include);

        EXPECT_EQ(r.n_atoms(), 2u);
        EXPECT_EQ(r.name, "LYS");
        EXPECT_EQ(r.seqnum, 9);

        EXPECT_EQ(r[0].name, " C  ");
        EXPECT_NEAR(r[0].x, 26.127, 1e-4);
        EXPECT_NEAR(r[0].y, 13.434, 1e-4);
        EXPECT_NEAR(r[0].z, 69.495, 1e-4);

        EXPECT_EQ(r[1].name, " O  ");
        EXPECT_NEAR(r[1].x, 27.16, 1e-4);
        EXPECT_NEAR(r[1].y, 14.073, 1e-4);
        EXPECT_NEAR(r[1].z, 69.678, 1e-4);
    }

    // Case 2: Valid residue with N and H atoms
    {
        std::vector<std::string> atoms_to_include;
        atoms_to_include.push_back(" N  ");
        atoms_to_include.push_back(" H  ");

        coot::minimol::residue r = coot::backrub::make_residue_include_only_gemmi(next_res, atoms_to_include);

        EXPECT_EQ(r.n_atoms(), 1u);
        EXPECT_EQ(r.name, "GLY");
        EXPECT_EQ(r.seqnum, 11);

        EXPECT_EQ(r[0].name, " N  ");
        EXPECT_NEAR(r[0].x, 25.736, 1e-4);
        EXPECT_NEAR(r[0].y, 13.244, 1e-4);
        EXPECT_NEAR(r[0].z, 66.003, 1e-4);
    }

    // Case 3: Single-atom filter (only CA)
    {
        std::vector<std::string> atoms_to_include;
        atoms_to_include.push_back(" CA ");

        coot::minimol::residue r = coot::backrub::make_residue_include_only_gemmi(this_res, atoms_to_include);

        EXPECT_EQ(r.n_atoms(), 1u);
        EXPECT_EQ(r.name, "ILE");
        EXPECT_EQ(r.seqnum, 10);

        EXPECT_EQ(r[0].name, " CA ");
        EXPECT_NEAR(r[0].x, 27.175, 1e-4);
        EXPECT_NEAR(r[0].y, 12.237, 1e-4);
        EXPECT_NEAR(r[0].z, 67.71, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}