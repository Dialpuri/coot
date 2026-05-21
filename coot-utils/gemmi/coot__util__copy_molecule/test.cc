#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Test 1: Normal molecule copy (protein with 2 chains, 268 residues)
TEST(CopyMoleculeTest, CopyMoleculeFull) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    int n_chains_original = 0;
    int n_residues_original = 0;
    if (!st.models.empty()) {
        n_chains_original = st.models[0].chains.size();
        for (const auto& chain : st.models[0].chains) {
            n_residues_original += chain.residues.size();
        }
    }

    EXPECT_EQ(n_chains_original, 2);
    EXPECT_EQ(n_residues_original, 268);

    gemmi::Structure st_copy = coot::util::copy_molecule_gemmi(st);

    int n_chains_copy = 0;
    int n_residues_copy = 0;
    if (!st_copy.models.empty()) {
        n_chains_copy = st_copy.models[0].chains.size();
        for (const auto& chain : st_copy.models[0].chains) {
            n_residues_copy += chain.residues.size();
        }
    }

    EXPECT_EQ(n_chains_copy, 2);
    EXPECT_EQ(n_residues_copy, 268);
    EXPECT_EQ(n_residues_original, n_residues_copy);
}

// Test 2: Copy then verify atom coordinates match (subset test)
TEST(CopyMoleculeTest, CopyMoleculeCoords) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Structure st_copy = coot::util::copy_molecule_gemmi(st);

    if (!st.models.empty() && !st_copy.models.empty()) {
        const auto& chain = st.models[0].chains[0];
        const auto& chain_copy = st_copy.models[0].chains[0];

        if (!chain.residues.empty() && !chain_copy.residues.empty()) {
            const auto& res = chain.residues[0];
            const auto& res_copy = chain_copy.residues[0];

            if (!res.atoms.empty() && !res_copy.atoms.empty()) {
                const auto& atom = res.atoms[0];
                const auto& atom_copy = res_copy.atoms[0];

                double dx = atom.pos.x - atom_copy.pos.x;
                double dy = atom.pos.y - atom_copy.pos.y;
                double dz = atom.pos.z - atom_copy.pos.z;
                bool coords_match = (std::fabs(dx) < 1e-6 && std::fabs(dy) < 1e-6 && std::fabs(dz) < 1e-6);

                EXPECT_TRUE(coords_match);
            }
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}