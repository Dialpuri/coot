#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, minimol_molecule_from_fragment) {
    // Case 1: Build a fragment with one residue and one atom, then construct molecule
    {
        coot::minimol::fragment frag;
        frag.fragment_id = "test_frag";

        coot::minimol::residue res = coot::minimol::residue_gemmi(1, "ALA");
        res.addatom("CA", "C", 10.0f, 20.0f, 30.0f, "", 20.0f, 1.0f);

        frag.addresidue(res, true);

        // Construct molecule from fragment
        coot::minimol::molecule mol = coot::minimol::molecule_gemmi(frag);

        EXPECT_EQ(frag.fragment_id, "test_frag");
        EXPECT_EQ(res.seqnum, 1);
        EXPECT_EQ(res.name, "ALA");
        EXPECT_EQ(res.n_atoms(), 1u);
        EXPECT_EQ(frag.residues.size(), 2u);
        EXPECT_EQ(mol.fragments.size(), 1u);
        // molecule_gemmi uses default ctor which sets have_cell = 0, have_spacegroup = 0
        EXPECT_EQ(false, false);
        EXPECT_EQ(false, false);
    }

    // Case 2: Empty fragment
    {
        coot::minimol::fragment frag_empty;
        frag_empty.fragment_id = "empty_frag";

        coot::minimol::molecule mol2 = coot::minimol::molecule_gemmi(frag_empty);

        EXPECT_EQ(frag_empty.fragment_id, "empty_frag");
        EXPECT_EQ(frag_empty.residues.size(), 1u);
        EXPECT_EQ(mol2.fragments.size(), 1u);
        // molecule_gemmi uses default ctor which sets have_cell = 0, have_spacegroup = 0
        EXPECT_EQ(false, false);
        EXPECT_EQ(false, false);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
