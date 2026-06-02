#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: find a residue by chain ID and sequence number in model 0
static gemmi::Residue* find_residue(gemmi::Structure& st, char chain_id, int seqnum) {
    gemmi::Model& model = st.models[0];
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name != std::string(1, chain_id)) continue;
        for (gemmi::Residue& res : chain.residues) {
            if (res.seqid.num.value == seqnum) return &res;
        }
    }
    return nullptr;
}

TEST(OracleTestGemmi, copy_and_delete_hydrogens) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    // Case 1: residue WITH hydrogens (MET 1 has 19 atoms including H)
    {
        gemmi::Residue* res = find_residue(st, 'A', 1);
        ASSERT_NE(res, nullptr);

        EXPECT_STREQ(res->name.c_str(), "MET");

        int original_atoms = static_cast<int>(res->atoms.size());
        EXPECT_EQ(original_atoms, 19);

        gemmi::Residue* cleaned = coot::util::copy_and_delete_hydrogens_gemmi(res);
        ASSERT_NE(cleaned, nullptr);

        int cleaned_atoms = static_cast<int>(cleaned->atoms.size());
        EXPECT_EQ(cleaned_atoms, 8);
        EXPECT_EQ(original_atoms - cleaned_atoms, 11);

        delete cleaned;
    }

    // Case 2: residue with fewer hydrogens (GLY 11 has 7 atoms)
    {
        gemmi::Residue* res = find_residue(st, 'A', 11);
        ASSERT_NE(res, nullptr);

        EXPECT_STREQ(res->name.c_str(), "GLY");

        int original_atoms = static_cast<int>(res->atoms.size());
        EXPECT_EQ(original_atoms, 7);

        gemmi::Residue* cleaned = coot::util::copy_and_delete_hydrogens_gemmi(res);
        ASSERT_NE(cleaned, nullptr);

        int cleaned_atoms = static_cast<int>(cleaned->atoms.size());
        EXPECT_EQ(cleaned_atoms, 4);
        EXPECT_EQ(original_atoms - cleaned_atoms, 3);

        delete cleaned;
    }

    // Case 3: residue with no hydrogens (ALA 140 has 6 atoms, no H)
    {
        gemmi::Residue* res = find_residue(st, 'A', 140);
        ASSERT_NE(res, nullptr);

        EXPECT_STREQ(res->name.c_str(), "ALA");

        int original_atoms = static_cast<int>(res->atoms.size());
        EXPECT_EQ(original_atoms, 6);

        gemmi::Residue* cleaned = coot::util::copy_and_delete_hydrogens_gemmi(res);
        ASSERT_NE(cleaned, nullptr);

        int cleaned_atoms = static_cast<int>(cleaned->atoms.size());
        EXPECT_EQ(cleaned_atoms, 4);
        EXPECT_EQ(original_atoms - cleaned_atoms, 2);

        delete cleaned;
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}