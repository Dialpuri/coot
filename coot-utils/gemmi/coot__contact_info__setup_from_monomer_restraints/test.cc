#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, contact_info) {
    // Load PDB using gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    ASSERT_EQ(st.models.size(), 1u);

    coot::protein_geometry geom;
    geom.init_standard();

    gemmi::Model& model = st.models[0];

    // Case 1: Single residue selection - //A/10/*
    {
        std::vector<gemmi::Atom*> sel_atoms;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name != "A") continue;
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value != 10) continue;
                for (gemmi::Atom& atom : res.atoms) {
                    sel_atoms.push_back(&atom);
                }
            }
        }

        EXPECT_EQ(sel_atoms.size(), 8u);

        std::vector<std::pair<int, int>> contacts;
        coot::setup_from_monomer_restraints_gemmi(
            model, sel_atoms, 0, &geom, contacts);

        EXPECT_EQ(contacts.size(), 14u);

        EXPECT_EQ(contacts[0].first, 0);
        EXPECT_EQ(contacts[0].second, 1);
        EXPECT_EQ(contacts[1].first, 0);
        EXPECT_EQ(contacts[1].second, 1);
        EXPECT_EQ(contacts[2].first, 1);
        EXPECT_EQ(contacts[2].second, 2);
        EXPECT_EQ(contacts[3].first, 1);
        EXPECT_EQ(contacts[3].second, 4);
        EXPECT_EQ(contacts[4].first, 1);
        EXPECT_EQ(contacts[4].second, 2);
    }

    // Case 2: Entire chain A - //A/*/*
    {
        std::vector<gemmi::Atom*> sel_atoms;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name != "A") continue;
            for (gemmi::Residue& res : chain.residues) {
                for (gemmi::Atom& atom : res.atoms) {
                    sel_atoms.push_back(&atom);
                }
            }
        }

        EXPECT_EQ(sel_atoms.size(), 2106u);

        std::vector<std::pair<int, int>> contacts;
        coot::setup_from_monomer_restraints_gemmi(
            model, sel_atoms, 0, &geom, contacts);

        EXPECT_EQ(contacts.size(), 3778u);

        EXPECT_EQ(contacts[0].first, 0);
        EXPECT_EQ(contacts[0].second, 1);
        EXPECT_EQ(contacts[1].first, 0);
        EXPECT_EQ(contacts[1].second, 1);
        EXPECT_EQ(contacts[2].first, 1);
        EXPECT_EQ(contacts[2].second, 2);
        EXPECT_EQ(contacts[3].first, 1);
        EXPECT_EQ(contacts[3].second, 4);
        EXPECT_EQ(contacts[4].first, 1);
        EXPECT_EQ(contacts[4].second, 2);
    }

    // Case 3: Selection //B/10/* → 0 atoms
    {
        std::vector<gemmi::Atom*> sel_atoms;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name != "B") continue;
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value != 10) continue;
                for (gemmi::Atom& atom : res.atoms) {
                    sel_atoms.push_back(&atom);
                }
            }
        }

        EXPECT_EQ(sel_atoms.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}