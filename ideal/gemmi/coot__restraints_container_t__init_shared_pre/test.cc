#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// --- Case 1: Call init_shared_pre_gemmi with a valid mol pointer ---
TEST(OracleTestGemmi, init_shared_pre_valid) {
  // Setup: load a PDB to get a gemmi::Structure*
  gemmi::Structure mol_structure = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

  // --- Case 1: Call init_shared_pre_gemmi with a valid mol pointer ---
  {
    // Construct a restraints_container_t via default constructor
    coot::restraints_container_t rc;

    // Record state BEFORE
    bool needs_reset_before = rc.needs_reset;
    double geman_before = rc.geman_mcclure_alpha;
    bool cryo_before = rc.cryo_em_mode;
    unsigned int n_calls_before = rc.n_times_called;

    EXPECT_FALSE(needs_reset_before);
    EXPECT_NEAR(geman_before, 0.0, 1e-4);
    EXPECT_FALSE(cryo_before);
    EXPECT_EQ(n_calls_before, 0u);

    // Call the function
    rc.init_shared_pre_gemmi(&mol_structure);

    // Record state AFTER
    EXPECT_FALSE(rc.needs_reset);
    EXPECT_NEAR(rc.geman_mcclure_alpha, 0.2, 1e-4);
    EXPECT_TRUE(rc.cryo_em_mode);
    EXPECT_TRUE(rc.mol == &mol_structure);
    EXPECT_EQ(rc.n_times_called, 0u);
    EXPECT_EQ(rc.dist_crit_for_bonded_pairs, 3);
    EXPECT_NEAR(rc.lennard_jones_epsilon, 1.0, 1e-4);
    EXPECT_NEAR(rc.log_cosh_target_distance_scale_factor, 3000.0, 1e-4);
    EXPECT_EQ(rc.verbose_geometry_reporting, 1);
    EXPECT_FALSE(rc.do_numerical_gradients_flag);
    EXPECT_FALSE(rc.have_oxt_flag);
    EXPECT_EQ(rc.n_small_cycles_accumulator, 0);
    EXPECT_FALSE(rc.use_proportional_editing);
    EXPECT_NEAR(rc.pull_restraint_neighbour_displacement_max_radius, 10.0, 1e-4);
    EXPECT_FALSE(rc.convert_plane_restraints_to_improper_dihedral_restraints_flag);
    EXPECT_EQ(rc.n_threads, 0);
  }

  // --- Case 2: Call init_shared_pre_gemmi with nullptr (tests guard behavior) ---
  {
    coot::restraints_container_t rc2;
    rc2.init_shared_pre_gemmi(nullptr);
    EXPECT_TRUE(rc2.mol == nullptr);
    EXPECT_FALSE(rc2.needs_reset);
    EXPECT_TRUE(rc2.cryo_em_mode);
  }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}