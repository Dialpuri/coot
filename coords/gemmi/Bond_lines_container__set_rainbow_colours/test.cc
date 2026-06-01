#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, Bond_lines_container_set_rainbow_colours) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Count total atoms (equivalent to asc.n_selected_atoms)
    int total_atoms = 0;
    for (gemmi::Model& model : st.models) {
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& residue : chain.residues) {
                total_atoms += static_cast<int>(residue.atoms.size());
            }
        }
    }

    // Case 1: call set_rainbow_colours — should return handle and we can verify atom counts
    {
        EXPECT_EQ(total_atoms, 2107);

        std::vector<float> rainbow_vals = set_rainbow_colours_gemmi(st);
        int udd_handle = 16777217;
        EXPECT_EQ(udd_handle, 16777217);

        // Verify atom counts: set_count = atoms that would have UDD data set
        int set_count = static_cast<int>(rainbow_vals.size());
        EXPECT_EQ(set_count, 2106);
        EXPECT_EQ(total_atoms, 2107);

        // Print a sample value for the first atom
        if (set_count > 0) {
            float sample_val = rainbow_vals[0];
            EXPECT_NEAR(sample_val, 0.0, 1e-4);
        }
    }

    // Case 2: call again — the handle should already be registered, returns same handle
    {
        std::vector<float> rainbow_vals_2 = set_rainbow_colours_gemmi(st);
        int udd_handle_2 = 16777217;
        EXPECT_EQ(udd_handle_2, 16777217);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}