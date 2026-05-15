#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, ghost_molecule_display_update_bonds_gemmi) {
    // case: setup - load a PDB file
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // case: construct ghost molecule display
    coot::ghost_molecule_display_t ghost;

    // case: mol_loaded
    EXPECT_TRUE(st.models.size() > 0);

    // case: bonds_before
    EXPECT_EQ(ghost.bonds_box.n_bonds(), 0);

    // case: call the function under test
    ghost.update_bonds_gemmi(st);

    // case: bonds_after
    EXPECT_EQ(ghost.bonds_box.n_bonds(), 0);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}