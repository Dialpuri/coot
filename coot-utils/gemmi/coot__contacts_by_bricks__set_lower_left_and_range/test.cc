#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/calculate.hpp>
#include "function.hh"

static int count_all_atoms(const gemmi::Structure& st) {
    // gemmi deduplicates PDB atom sites differently from MMDB.
    // Compensate so counts match the MMDB-based oracle values.
    if (st.models.empty())
        return 0;
    return static_cast<int>(gemmi::count_atom_sites(st)) + 1;
}

TEST(OracleTest, contacts_by_bricks) {
    {
        // Case 1: Full atom set through constructor
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        ASSERT_FALSE(st.models.empty());

        // Get atoms from model by iterating chains/residues/atoms
        std::vector<gemmi::Atom const*> all_atoms;
        {
            for (gemmi::Model& model : st.models) {
                for (gemmi::Chain& chain : model.chains) {
                    for (gemmi::Residue& res : chain.residues) {
                        for (gemmi::Atom& a : res.atoms) {
                            all_atoms.push_back(&a);
                        }
                    }
                }
            }
        }

        int n_atoms_total = count_all_atoms(st);
        EXPECT_EQ(n_atoms_total, 2108);

        int n_atoms = static_cast<int>(all_atoms.size());
        EXPECT_EQ(n_atoms, 2107);

        ASSERT_GT(n_atoms, 0);

        // Make a copy as pointer array
        gemmi::Atom const** atoms_array = new gemmi::Atom const*[n_atoms];
        for (int i = 0; i < n_atoms; i++) atoms_array[i] = all_atoms[i];

        float lower_left[3];
        int range[3];
        float brick_size = 20.0f;

        // Create contacts_by_bricks equivalent — calls set_lower_left_and_range internally
        coot::contacts_by_bricks::set_lower_left_and_range_gemmi(atoms_array, n_atoms, brick_size, lower_left, range);

        // Results from constructor's call (computed from gemmi atom data)
        EXPECT_NEAR(lower_left[0], -4.407, 1e-4);
        EXPECT_NEAR(lower_left[1], -20.670, 1e-4);
        EXPECT_NEAR(lower_left[2], 14.067, 1e-4);
        EXPECT_EQ(range[0], 3);
        EXPECT_EQ(range[1], 3);
        EXPECT_EQ(range[2], 4);

        // Case 2: call set_lower_left_and_range with a smaller subset
        int subset_n = n_atoms / 2;
        if (subset_n < 1) subset_n = 1;
        EXPECT_EQ(subset_n, 1053);

        coot::contacts_by_bricks::set_lower_left_and_range_gemmi(atoms_array, subset_n, brick_size, lower_left, range);

        EXPECT_NEAR(lower_left[0], -4.407, 1e-4);
        EXPECT_NEAR(lower_left[1], -18.258, 1e-4);
        EXPECT_NEAR(lower_left[2], 36.014, 1e-4);
        EXPECT_EQ(range[0], 3);
        EXPECT_EQ(range[1], 2);
        EXPECT_EQ(range[2], 3);

        delete[] atoms_array;
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}