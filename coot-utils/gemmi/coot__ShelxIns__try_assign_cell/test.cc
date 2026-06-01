#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include <cmath>

TEST(OracleTest, try_assign_cell) {
    // Load the PDB into a gemmi::Structure
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Case 1: Fresh ShelxIns with have_cell_flag = false, try to assign cell
    {
        coot::ShelxIns shelx;
        shelx.have_cell_flag = false;

        bool result = shelx.try_assign_cell_gemmi(st);

        EXPECT_TRUE(result);
        EXPECT_TRUE(shelx.have_cell_flag);
        EXPECT_NEAR(shelx.cell.a(), 53.512, 1e-3);
        EXPECT_NEAR(shelx.cell.b(), 71.311, 1e-3);
        EXPECT_NEAR(shelx.cell.c(), 72.101, 1e-3);
        EXPECT_NEAR(shelx.cell.alpha(), M_PI / 2.0, 1e-6);
        EXPECT_NEAR(shelx.cell.beta(), M_PI / 2.0, 1e-6);
        EXPECT_NEAR(shelx.cell.gamma(), M_PI / 2.0, 1e-6);
    }

    // Case 2: Pre-set have_cell_flag = true — should return true immediately
    {
        coot::ShelxIns shelx;
        shelx.have_cell_flag = true;

        bool result = shelx.try_assign_cell_gemmi(st);

        EXPECT_TRUE(result);
        EXPECT_TRUE(shelx.have_cell_flag);
    }

    // Case 3: Omitted — original tested crash on nullptr mmdb::Manager*.
    // With const gemmi::Structure& (reference), nullptr is not representable,
    // so this DEATH test does not translate.
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}