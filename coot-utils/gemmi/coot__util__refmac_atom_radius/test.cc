#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <cmath>
#include "function.hh"

TEST(OracleTest, refmac_atom_radius) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    ASSERT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];

    // Get first chain
    ASSERT_FALSE(model.chains.empty());
    gemmi::Chain& chain = model.chains[0];

    // Get first residue
    ASSERT_FALSE(chain.residues.empty());
    gemmi::Residue& res = chain.residues[0];

    // Case 1: first atom of first residue
    {
        ASSERT_FALSE(res.atoms.empty());
        gemmi::Atom& at = res.atoms[0];
        std::string name = at.name;
        double bf = at.b_iso;
        std::string resname = res.name;
        std::string chainid = chain.name;
        int seqnum = res.seqid.num.value;

        double radius = coot::util::refmac_atom_radius_gemmi(at);

        EXPECT_EQ(chainid, "A");
        EXPECT_EQ(seqnum, 1);
        EXPECT_EQ(resname, "MET");
        EXPECT_NEAR(bf, 65.05, 1e-4);
        EXPECT_NEAR(radius, 2.99849, 1e-4);
    }

    // Case 2: another residue further down the chain
    {
        gemmi::Residue* res2 = nullptr;
        if (chain.residues.size() > 10) {
            res2 = &chain.residues[10];
        } else if (!chain.residues.empty()) {
            res2 = &chain.residues[1];
        }
        ASSERT_NE(res2, nullptr);

        gemmi::Atom* at = nullptr;
        if (!res2->atoms.empty()) {
            at = &res2->atoms[0];
        } else {
            at = &res.atoms[0];
        }
        ASSERT_NE(at, nullptr);

        std::string name = at->name;
        double bf = at->b_iso;
        std::string resname = res2->name;
        std::string chainid = chain.name;
        int seqnum = res2->seqid.num.value;

        double radius = coot::util::refmac_atom_radius_gemmi(*at);

        EXPECT_EQ(chainid, "A");
        EXPECT_EQ(seqnum, 11);
        EXPECT_EQ(resname, "GLY");
        EXPECT_NEAR(bf, 46.5, 1e-4);
        EXPECT_NEAR(radius, 2.66673, 1e-4);
    }

    // Case 3: atom with high B-factor (last residue in chain)
    {
        ASSERT_FALSE(chain.residues.empty());
        gemmi::Residue& res3 = chain.residues[chain.residues.size() - 1];
        gemmi::Atom* at = nullptr;
        if (!res3.atoms.empty()) {
            at = &res3.atoms[0];
        } else {
            at = &res.atoms[0];
        }
        ASSERT_NE(at, nullptr);

        std::string name = at->name;
        double bf = at->b_iso;
        std::string resname = res3.name;
        std::string chainid = chain.name;
        int seqnum = res3.seqid.num.value;

        double radius = coot::util::refmac_atom_radius_gemmi(*at);

        EXPECT_EQ(chainid, "A");
        EXPECT_EQ(seqnum, 298);
        EXPECT_EQ(resname, "LEU");
        EXPECT_NEAR(bf, 42.59, 1e-4);
        EXPECT_NEAR(radius, 2.59138, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}