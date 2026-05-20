#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, MoleculeCell) {
    // Create an empty molecule
    coot::minimol::molecule mol;

    // Print initial state
    std::vector<float> initial_cell = mol.get_cell();
    EXPECT_EQ(initial_cell.size(), 0u);
    EXPECT_EQ(mol.have_cell, 0);

    // Test case 1: Set cell with float array (standard case)
    {
        float a[6] = {10.0, 20.0, 30.0, 90.0, 90.0, 90.0};
        mol.set_cell(a);

        std::vector<float> new_cell = mol.get_cell();
        EXPECT_EQ(new_cell.size(), 6u);
        EXPECT_EQ(mol.have_cell, 1);
        EXPECT_EQ(new_cell[0], 10.0f);
        EXPECT_EQ(new_cell[1], 20.0f);
        EXPECT_EQ(new_cell[2], 30.0f);
        EXPECT_EQ(new_cell[3], 90.0f);
        EXPECT_EQ(new_cell[4], 90.0f);
        EXPECT_EQ(new_cell[5], 90.0f);
    }

    // Test case 2: Set cell again with different values (mutation check)
    {
        float b[6] = {50.0, 60.0, 70.0, 120.0, 90.0, 90.0};
        mol.set_cell(b);

        std::vector<float> updated_cell = mol.get_cell();
        EXPECT_EQ(updated_cell.size(), 6u);
        EXPECT_EQ(updated_cell[0], 50.0f);
        EXPECT_EQ(updated_cell[1], 60.0f);
        EXPECT_EQ(updated_cell[2], 70.0f);
        EXPECT_EQ(updated_cell[3], 120.0f);
        EXPECT_EQ(updated_cell[4], 90.0f);
        EXPECT_EQ(updated_cell[5], 90.0f);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}