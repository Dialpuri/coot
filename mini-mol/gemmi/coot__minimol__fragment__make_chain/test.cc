#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, fragment_make_chain) {
    // Case 1: fragment with residues containing atoms
    {
        coot::minimol::fragment frag("A");

        // Add residue 10 with a CA atom
        coot::minimol::residue res10(10, "ALA");
        res10.addatom("CA", " C", 0.0f, 0.0f, 0.0f, "", 20.0f, 1.0f);
        res10.addatom("N",  " N", 1.0f, 0.0f, 0.0f, "", 20.0f, 1.0f);
        frag.addresidue(res10, true);

        // Add residue 11 with a CA atom
        coot::minimol::residue res11(11, "GLY");
        res11.addatom("CA", " C", 2.0f, 0.0f, 0.0f, "", 25.0f, 1.0f);
        frag.addresidue(res11, true);

        EXPECT_EQ(frag.fragment_id, "A");
        EXPECT_EQ(frag.min_res_no(), 1);
        EXPECT_EQ(frag.max_residue_number(), 11);

        auto chain_p = make_chain_gemmi(frag);

        EXPECT_NE(chain_p, nullptr);
        EXPECT_EQ(chain_p->name, "A");
        EXPECT_EQ(chain_p->residues.size(), 2u);

        if (chain_p) {
            // residue_0
            const gemmi::Residue& r0 = chain_p->residues[0];
            EXPECT_EQ(chain_p->name, "A");
            EXPECT_EQ(r0.seqid.num.value, 10);
            EXPECT_EQ(r0.name, "ALA");
            EXPECT_EQ(r0.atoms.size(), 2u);

            // residue_1
            const gemmi::Residue& r1 = chain_p->residues[1];
            EXPECT_EQ(chain_p->name, "A");
            EXPECT_EQ(r1.seqid.num.value, 11);
            EXPECT_EQ(r1.name, "GLY");
            EXPECT_EQ(r1.atoms.size(), 1u);
        }
    }

    // Case 2: default-constructed empty fragment (edge case — no real residues)
    {
        coot::minimol::fragment frag_empty("B");

        EXPECT_EQ(frag_empty.fragment_id, "B");
        EXPECT_EQ(frag_empty.min_res_no(), 1);
        EXPECT_EQ(frag_empty.max_residue_number(), 0);

        auto chain_p = make_chain_gemmi(frag_empty);

        EXPECT_NE(chain_p, nullptr);
        EXPECT_EQ(chain_p->name, "B");
        EXPECT_EQ(chain_p->residues.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}