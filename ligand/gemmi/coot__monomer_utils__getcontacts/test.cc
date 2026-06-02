#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, getcontacts) {
    // Load PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    EXPECT_EQ(st.models.size(), 1u);

    gemmi::Model& model = st.models[0];

    // Case 1: Select atoms from chain A (should find contacts within chain A)
    {
        std::vector<gemmi::Atom*> selected_atoms;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue& res : chain.residues) {
                    for (gemmi::Atom& atom : res.atoms) {
                        selected_atoms.push_back(&atom);
                    }
                }
            }
        }

        EXPECT_EQ(selected_atoms.size(), 2106u);

        std::vector<std::pair<int, int>> contacts = coot::monomer_utils::getcontacts_gemmi(model, selected_atoms);

        EXPECT_EQ(contacts.size(), 4304u);

        // First few contacts
        EXPECT_GE(contacts.size(), 5u);
        {
            auto cp = contacts[0];
            EXPECT_EQ(cp.first, 0);
            EXPECT_EQ(cp.second, 1);
        }
        {
            auto cp = contacts[1];
            EXPECT_EQ(cp.first, 1);
            EXPECT_EQ(cp.second, 0);
        }
        {
            auto cp = contacts[2];
            EXPECT_EQ(cp.first, 1);
            EXPECT_EQ(cp.second, 4);
        }
        {
            auto cp = contacts[3];
            EXPECT_EQ(cp.first, 1);
            EXPECT_EQ(cp.second, 2);
        }
        {
            auto cp = contacts[4];
            EXPECT_EQ(cp.first, 2);
            EXPECT_EQ(cp.second, 1);
        }
    }

    // Case 2: Select atoms from chain B (should find contacts within chain B)
    {
        std::vector<gemmi::Atom*> selected_atoms;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "B") {
                for (gemmi::Residue& res : chain.residues) {
                    for (gemmi::Atom& atom : res.atoms) {
                        selected_atoms.push_back(&atom);
                    }
                }
            }
        }

        EXPECT_EQ(selected_atoms.size(), 1u);

        std::vector<std::pair<int, int>> contacts = coot::monomer_utils::getcontacts_gemmi(model, selected_atoms);

        EXPECT_EQ(contacts.size(), 0u);
    }

    // Case 3: Select a single residue (small selection, edge case)
    // Oracle output shows only case3_n_selected_atoms: 0 and no n_contacts output,
    // indicating getcontacts crashes with 0 selected atoms.
    // MMDB: STYPE_RESIDUE selects a RESIDUE (not ATOMS), so atom_selection_container_t
    // reports n_selected_atoms=0. We replicate this by creating an empty selection.
    {
        std::vector<gemmi::Atom*> selected_atoms; // empty — matches MMDB STYPE_RESIDUE behavior

        EXPECT_EQ(selected_atoms.size(), 0u);

        // getcontacts should handle 0 selected atoms without crashing
        std::vector<std::pair<int, int>> contacts = coot::monomer_utils::getcontacts_gemmi(model, selected_atoms);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}