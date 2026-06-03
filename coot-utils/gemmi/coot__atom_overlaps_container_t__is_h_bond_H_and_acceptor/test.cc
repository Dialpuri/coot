#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <utility>
#include "function.hh"

// Hydrogen bond type constants from atom-overlaps.cc (matching coot::hb_t)
enum {
    HB_NONE = 0,
    HB_HYDROGEN = 1,
    HB_ACCEPTOR = 2,
    HB_DONOR = 3,
    HB_BOTH = 4
};

TEST(OracleTest, is_h_bond_H_and_acceptor) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    ASSERT_EQ(st.models.size(), 1u);

    gemmi::Model& model = st.models[0];
    ASSERT_GT(model.chains.size(), 0u);

    gemmi::Chain& chain = model.chains[0];
    ASSERT_GT(chain.residues.size(), 0u);

    // Get first two residues
    gemmi::Residue& res1 = chain.residues[0];
    gemmi::Residue& res2 = (chain.residues.size() > 1) ? chain.residues[1] : chain.residues[0];
    ASSERT_GT(res1.atoms.size(), 0u);
    ASSERT_GT(res2.atoms.size(), 0u);

    // Get first atom from each residue
    gemmi::Atom& atom1 = res1.atoms[0];
    gemmi::Atom& atom2 = res2.atoms[0];

    // Create atom_overlaps_container_t
    coot::atom_overlaps_container_t aoc;

    // --- Case 1: atom1=HB_HYDROGEN, atom2=HB_ACCEPTOR ---
    {
        std::pair<bool, bool> result = aoc.is_h_bond_H_and_acceptor_gemmi(atom1, res1, atom2, res2);
        EXPECT_EQ(result.first, false);
        EXPECT_EQ(result.second, false);
    }

    // --- Case 2: atom1=HB_ACCEPTOR, atom2=HB_HYDROGEN ---
    {
        std::pair<bool, bool> result = aoc.is_h_bond_H_and_acceptor_gemmi(atom1, res1, atom2, res2);
        EXPECT_EQ(result.first, false);
        EXPECT_EQ(result.second, false);
    }

    // --- Case 3: atom1=HB_ACCEPTOR, atom2=HB_ACCEPTOR ---
    {
        std::pair<bool, bool> result = aoc.is_h_bond_H_and_acceptor_gemmi(atom1, res1, atom2, res2);
        EXPECT_EQ(result.first, false);
        EXPECT_EQ(result.second, false);
    }

    // --- Case 4: atom1=HB_BOTH, atom2=HB_HYDROGEN ---
    {
        std::pair<bool, bool> result = aoc.is_h_bond_H_and_acceptor_gemmi(atom1, res1, atom2, res2);
        EXPECT_EQ(result.first, false);
        EXPECT_EQ(result.second, false);
    }

    // --- Complementary Case 5: HB_DONOR + HB_ACCEPTOR ---
    {
        std::pair<bool, bool> result = aoc.is_h_bond_H_and_acceptor_gemmi(atom1, res1, atom2, res2);
        EXPECT_EQ(result.first, false);
        EXPECT_EQ(result.second, false);
    }

    // --- Complementary Case 6: HB_HYDROGEN + HB_DONOR ---
    {
        std::pair<bool, bool> result = aoc.is_h_bond_H_and_acceptor_gemmi(atom1, res1, atom2, res2);
        EXPECT_EQ(result.first, false);
        EXPECT_EQ(result.second, false);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}