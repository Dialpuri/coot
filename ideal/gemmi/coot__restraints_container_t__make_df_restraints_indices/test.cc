#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, make_df_restraints_indices_gemmi) {
    {
        coot::restraints_container_t_gemmi restraints;

        EXPECT_EQ(restraints.size(), 0u);
        EXPECT_EQ(restraints.get_n_atoms(), 0u);

        // n_threads defaults to 0, function should handle gracefully (fallback to 1)
        restraints.n_threads = 0;

        restraints.make_df_restraints_indices_gemmi();

        // With n_threads==0, n_r_s becomes 1
        EXPECT_EQ(restraints.restraints_indices.size(), 1u);
        EXPECT_EQ(restraints.df_by_thread_results.size(), 1u);
        EXPECT_EQ(restraints.df_by_thread_atom_indices.size(), 1u);

        // Empty restraints -> empty indices
        EXPECT_EQ(restraints.restraints_indices[0].size(), 0u);
        EXPECT_EQ(restraints.df_by_thread_atom_indices[0].size(), 0u);
    }

    // Test with n_threads == 2
    {
        coot::restraints_container_t_gemmi restraints;
        restraints.n_threads = 2;

        restraints.make_df_restraints_indices_gemmi();

        EXPECT_EQ(restraints.restraints_indices.size(), 2u);
        EXPECT_EQ(restraints.df_by_thread_results.size(), 2u);
        EXPECT_EQ(restraints.df_by_thread_atom_indices.size(), 2u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
