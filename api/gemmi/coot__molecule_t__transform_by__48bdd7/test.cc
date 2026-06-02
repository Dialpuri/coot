#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, transform_by) {
    // Load PDB file into gemmi structure
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Get first atom coordinates BEFORE transform
    gemmi::Model &model = st.models[0];
    gemmi::Chain &chain = model.chains[0];
    gemmi::Residue &res = chain.residues[0];
    gemmi::Atom &first_atom = res.atoms[0];

    double x_before = first_atom.pos.x;
    double y_before = first_atom.pos.y;
    double z_before = first_atom.pos.z;

    EXPECT_NEAR(x_before, 12.334, 1e-4);
    EXPECT_NEAR(y_before, 2.772, 1e-4);
    EXPECT_NEAR(z_before, 74.951, 1e-4);

    // Case 1: Identity matrix (no change)
    std::array<std::array<double, 4>, 4> identity;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            identity[i][j] = (i == j) ? 1.0 : 0.0;

    coot::transform_by_gemmi(st, identity);

    double x_after_identity = first_atom.pos.x;
    double y_after_identity = first_atom.pos.y;
    double z_after_identity = first_atom.pos.z;

    EXPECT_NEAR(x_after_identity, 12.334, 1e-4);
    EXPECT_NEAR(y_after_identity, 2.772, 1e-4);
    EXPECT_NEAR(z_after_identity, 74.951, 1e-4);

    // Case 2: Translation by 10.0 Angstroms in X
    std::array<std::array<double, 4>, 4> translate_x;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            translate_x[i][j] = (i == j) ? 1.0 : 0.0;
    translate_x[0][3] = 10.0;

    coot::transform_by_gemmi(st, translate_x);

    double x_after_translate = first_atom.pos.x;
    double y_after_translate = first_atom.pos.y;
    double z_after_translate = first_atom.pos.z;

    EXPECT_NEAR(x_after_translate, 22.334, 1e-4);
    EXPECT_NEAR(y_after_translate, 2.772, 1e-4);
    EXPECT_NEAR(z_after_translate, 74.951, 1e-4);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}