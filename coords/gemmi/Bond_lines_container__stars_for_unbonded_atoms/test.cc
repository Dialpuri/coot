#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <unordered_set>
#include "function.hh"

TEST(OracleTest, stars_for_unbonded_atoms) {
    // Case 1: Empty bonds[2] — addBond guard rejects new bonds (same as original)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        gemmi::Model& model = st.models[0];

        // Build set of all atoms (all "unbonded")
        std::unordered_set<const gemmi::Atom*> unbonded;
        for (auto& chain : model.chains)
            for (auto& res : chain.residues)
                for (auto& atom : res.atoms)
                    unbonded.insert(&atom);

        // Count total atoms — gemmi reads 2107 (MMDB read 2108 including HETATM water)
        size_t total_atoms = 0;
        for (auto& chain : model.chains)
            for (auto& res : chain.residues)
                total_atoms += res.atoms.size();
        EXPECT_EQ(total_atoms, 2107);

        // Create Bond_lines_container_gemmi with indices 0,1 populated
        Bond_lines_container_gemmi bonds;
        Bond_lines a;
        bonds.bonds.push_back(a); // index 0
        bonds.bonds.push_back(a); // index 1

        int bonds_before = bonds.bonds[1].size();
        EXPECT_EQ(bonds_before, 0);

        EXPECT_NO_THROW(stars_for_unbonded_atoms_gemmi(bonds, model, unbonded));

        int bonds_after = bonds.bonds[1].size();
        EXPECT_EQ(bonds_after, 0);
        EXPECT_EQ(bonds_after - bonds_before, 0);
    }

    // Case 2: Pre-populate bonds[2] so addBond guard passes — core logic exercised
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        gemmi::Model& model = st.models[0];

        // Build set of all atoms (all "unbonded")
        std::unordered_set<const gemmi::Atom*> unbonded;
        for (auto& chain : model.chains)
            for (auto& res : chain.residues)
                for (auto& atom : res.atoms)
                    unbonded.insert(&atom);

        // Count total atoms
        int total_atoms = 0;
        for (auto& chain : model.chains)
            for (auto& res : chain.residues)
                total_atoms += static_cast<int>(res.atoms.size());
        EXPECT_GT(total_atoms, 0);

        // Start with empty container but ensure indices 0, 1, 2 exist
        Bond_lines_container_gemmi bonds;
        Bond_lines a;
        bonds.bonds.push_back(a);  // index 0
        bonds.bonds.push_back(a);  // index 1
        bonds.bonds.push_back(a);  // index 2

        // Pre-populate bonds[2] with one dummy entry so addBond guard passes
        bonds.bonds[2].add_bond(coot::CartesianPair(coot::Cartesian(0., 0., 0.), coot::Cartesian(0., 0., 0.)),
                                graphics_line_t::SINGLE, true, true, -1, -1, -1);

        // Count total bonds before
        int total_bonds_before = 0;
        for (const auto &bl : bonds.bonds) total_bonds_before += bl.size();

        // Call the target function
        EXPECT_NO_THROW(stars_for_unbonded_atoms_gemmi(bonds, model, unbonded));

        // Count total bonds after
        int total_bonds_after = 0;
        for (const auto &bl : bonds.bonds) total_bonds_after += bl.size();

        // Each NO_BOND atom gets 3 star bonds (x, y, z axes)
        int expected_new = 3 * total_atoms;
        int actual_new = total_bonds_after - total_bonds_before;
        EXPECT_EQ(actual_new, expected_new);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}