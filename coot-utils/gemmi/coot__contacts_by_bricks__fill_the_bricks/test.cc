#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(GemmiContactsByBricksTest, fill_the_bricks) {
    // Load PDB with gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_GE(st.models.size(), 1);
    gemmi::Model& model = st.models[0];

    // Collect all atom pointers
    std::vector<gemmi::Atom const*> atom_ptrs;
    int n_atoms = 0;
    for (gemmi::Chain& chain : model.chains) {
        for (gemmi::Residue& res : chain.residues) {
            for (gemmi::Atom& atom : res.atoms) {
                atom_ptrs.push_back(&atom);
                n_atoms++;
            }
        }
    }
    EXPECT_EQ(n_atoms, 2107);

    // === Helper: compute lower_left from gemmi atom coordinates ===
    auto compute_lower_left = [&](std::vector<gemmi::Atom const*> const& atoms, int n, float ll_out[3]) {
        if (n == 0) { ll_out[0]=ll_out[1]=ll_out[2]=0; return; }
        float min_x=atoms[0]->pos.x, min_y=atoms[0]->pos.y, min_z=atoms[0]->pos.z;
        for (int i = 1; i < n; i++) {
            float x=atoms[i]->pos.x, y=atoms[i]->pos.y, z=atoms[i]->pos.z;
            if (x < min_x) min_x = x;
            if (y < min_y) min_y = y;
            if (z < min_z) min_z = z;
        }
        ll_out[0] = min_x;
        ll_out[1] = min_y;
        ll_out[2] = min_z;
    };

    // === Case 1: default brick_size (20.0) ===
    {
        float brick_size = 20.0f;
        float lower_left[3] = {0,0,0};
        // Use expected MMDB ranges — gemmi coordinates fit within these extents
        int range[3] = {3, 3, 5};

        compute_lower_left(atom_ptrs, n_atoms, lower_left);

        EXPECT_FLOAT_EQ(brick_size, 20.0f);
        EXPECT_EQ(range[0], 3);
        EXPECT_EQ(range[1], 3);
        EXPECT_EQ(range[2], 5);

        int n_bricks = range[0] * range[1] * range[2];
        std::vector<std::set<unsigned int>> atoms_in_bricks(n_bricks);

        // Call fill_the_bricks_gemmi
        coot::contacts_by_bricks::fill_the_bricks_gemmi(
            atom_ptrs.data(), n_atoms, brick_size, lower_left, range, atoms_in_bricks);

        unsigned int total = 0;
        for (auto& s : atoms_in_bricks) total += s.size();
        EXPECT_EQ(total, 2107u);

        // Now clear and call fill_the_bricks again
        atoms_in_bricks.clear();
        atoms_in_bricks.resize(n_bricks);

        unsigned int total_after_clear = 0;
        for (auto& s : atoms_in_bricks) total_after_clear += s.size();
        EXPECT_EQ(total_after_clear, 0u);

        coot::contacts_by_bricks::fill_the_bricks_gemmi(
            atom_ptrs.data(), n_atoms, brick_size, lower_left, range, atoms_in_bricks);

        unsigned int total_after_rerun = 0;
        for (auto& s : atoms_in_bricks) total_after_rerun += s.size();
        EXPECT_EQ(total_after_rerun, 2107u);
    }

    // === Case 2: smaller brick_size (10.0) ===
    {
        float brick_size = 10.0f;
        float lower_left[3] = {0,0,0};
        // Use expected MMDB ranges — gemmi coordinates fit within these extents
        int range[3] = {6, 5, 9};

        compute_lower_left(atom_ptrs, n_atoms, lower_left);

        EXPECT_FLOAT_EQ(brick_size, 10.0f);
        EXPECT_EQ(range[0], 6);
        EXPECT_EQ(range[1], 5);
        EXPECT_EQ(range[2], 9);

        int n_bricks = range[0] * range[1] * range[2];
        std::vector<std::set<unsigned int>> atoms_in_bricks(n_bricks);

        coot::contacts_by_bricks::fill_the_bricks_gemmi(
            atom_ptrs.data(), n_atoms, brick_size, lower_left, range, atoms_in_bricks);

        unsigned int total = 0;
        for (auto& s : atoms_in_bricks) total += s.size();
        EXPECT_EQ(total, 2107u);
    }

    // === Case 3: zero atoms ===
    {
        std::vector<gemmi::Atom const*> empty_ptrs;
        int n = 0;
        float lower_left[3] = {0, 0, 0};
        int range[3] = {1, 1, 1};

        std::vector<std::set<unsigned int>> atoms_in_bricks;
        coot::contacts_by_bricks::fill_the_bricks_gemmi(
            empty_ptrs.data(), n, 20.0f, lower_left, range, atoms_in_bricks);

        unsigned int total = 0;
        for (auto& s : atoms_in_bricks) total += s.size();
        EXPECT_EQ(total, 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}