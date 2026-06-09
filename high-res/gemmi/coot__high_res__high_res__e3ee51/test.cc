#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, high_res) {
    // Case 1: valid protein molecule
    {
        coot::minimol::molecule mol;
        int load_result = mol.read_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        EXPECT_EQ(mol.get_number_of_atoms(), 2107);

        clipper::Coord_orth centre = mol.centre();

        EXPECT_NEAR(centre.x(), 29.4408, 1e-4);
        EXPECT_NEAR(centre.y(), 5.04435, 1e-4);
        EXPECT_NEAR(centre.z(), 50.2402, 1e-4);

        gemmi::Structure structure = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        coot::minimol::molecule globular_mol = coot::high_res::high_res_gemmi(mol, centre, structure);

        EXPECT_EQ(globular_mol.get_number_of_atoms(), 2107);
    }

    // Case 2: empty molecule (edge case)
    {
        coot::minimol::molecule empty_mol;
        clipper::Coord_orth zero_centre(0.0, 0.0, 0.0);

        EXPECT_EQ(empty_mol.get_number_of_atoms(), 0);
        EXPECT_NEAR(zero_centre.x(), 0.0, 1e-4);
        EXPECT_NEAR(zero_centre.y(), 0.0, 1e-4);
        EXPECT_NEAR(zero_centre.z(), 0.0, 1e-4);

        gemmi::Structure structure;

        coot::minimol::molecule globular_mol = coot::high_res::high_res_gemmi(empty_mol, zero_centre, structure);

        EXPECT_EQ(globular_mol.get_number_of_atoms(), 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
