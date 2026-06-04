#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, molecule_extents_t_gemmi) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty()) << "Failed to load PDB";
    gemmi::Model& model = st.models[0];

    // Count atoms for verification (equivalent to atom_sel.n_selected_atoms)
    int total_atoms = 0;
    for (const gemmi::Chain& chain : model.chains)
        for (const gemmi::Residue& res : chain.residues)
            total_atoms += static_cast<int>(res.atoms.size());

    // Case 1: Full molecule atom_sel with expansion size 20.0
    {
        float expansion_size = 20.0f;

        EXPECT_EQ(total_atoms, 2107);

        molecule_extents_t extents = molecule_extents_t::molecule_extents_t_gemmi(model, expansion_size);

        EXPECT_NEAR(extents.centre.x(), 26.481f, 1e-4);
        EXPECT_NEAR(extents.centre.y(), 5.563f, 1e-4);
        EXPECT_NEAR(extents.centre.z(), 51.505f, 1e-4);

        EXPECT_NEAR(extents.left.x(), -18.407f, 1e-4);
        EXPECT_NEAR(extents.left.y(), 5.563f, 1e-4);
        EXPECT_NEAR(extents.left.z(), 51.505f, 1e-4);

        EXPECT_NEAR(extents.right.x(), 71.369f, 1e-4);
        EXPECT_NEAR(extents.right.y(), 5.563f, 1e-4);
        EXPECT_NEAR(extents.right.z(), 51.505f, 1e-4);

        EXPECT_NEAR(extents.front.x(), 26.481f, 1e-4);
        EXPECT_NEAR(extents.front.y(), 5.563f, 1e-4);
        EXPECT_NEAR(extents.front.z(), 0.0669994f, 1e-4);

        EXPECT_NEAR(extents.back.x(), 26.481f, 1e-4);
        EXPECT_NEAR(extents.back.y(), 5.563f, 1e-4);
        EXPECT_NEAR(extents.back.z(), 102.943f, 1e-4);

        EXPECT_NEAR(extents.bottom.x(), 26.481f, 1e-4);
        EXPECT_NEAR(extents.bottom.y(), -34.67f, 1e-4);
        EXPECT_NEAR(extents.bottom.z(), 51.505f, 1e-4);

        EXPECT_NEAR(extents.top.x(), 26.481f, 1e-4);
        EXPECT_NEAR(extents.top.y(), 45.796f, 1e-4);
        EXPECT_NEAR(extents.top.z(), 51.505f, 1e-4);

        EXPECT_NEAR(extents.expansion_size_, 20.0f, 1e-4);
    }

    // Case 2: Same selection with expansion size 1.0
    {
        float expansion_size = 1.0f;

        EXPECT_EQ(total_atoms, 2107);

        molecule_extents_t extents = molecule_extents_t::molecule_extents_t_gemmi(model, expansion_size);

        EXPECT_NEAR(extents.centre.x(), 26.481f, 1e-4);
        EXPECT_NEAR(extents.centre.y(), 5.563f, 1e-4);
        EXPECT_NEAR(extents.centre.z(), 51.505f, 1e-4);

        EXPECT_NEAR(extents.left.x(), 0.593f, 1e-4);
        EXPECT_NEAR(extents.left.y(), 5.563f, 1e-4);
        EXPECT_NEAR(extents.left.z(), 51.505f, 1e-4);

        EXPECT_NEAR(extents.right.x(), 52.369f, 1e-4);
        EXPECT_NEAR(extents.right.y(), 5.563f, 1e-4);
        EXPECT_NEAR(extents.right.z(), 51.505f, 1e-4);

        EXPECT_NEAR(extents.front.x(), 26.481f, 1e-4);
        EXPECT_NEAR(extents.front.y(), 5.563f, 1e-4);
        EXPECT_NEAR(extents.front.z(), 19.067f, 1e-4);

        EXPECT_NEAR(extents.back.x(), 26.481f, 1e-4);
        EXPECT_NEAR(extents.back.y(), 5.563f, 1e-4);
        EXPECT_NEAR(extents.back.z(), 83.943f, 1e-4);

        EXPECT_NEAR(extents.bottom.x(), 26.481f, 1e-4);
        EXPECT_NEAR(extents.bottom.y(), -15.67f, 1e-4);
        EXPECT_NEAR(extents.bottom.z(), 51.505f, 1e-4);

        EXPECT_NEAR(extents.top.x(), 26.481f, 1e-4);
        EXPECT_NEAR(extents.top.y(), 26.796f, 1e-4);
        EXPECT_NEAR(extents.top.z(), 51.505f, 1e-4);

        EXPECT_NEAR(extents.expansion_size_, 1.0f, 1e-4);
    }

    // Case 3: Zero expansion size
    {
        float expansion_size = 0.0f;

        EXPECT_EQ(total_atoms, 2107);

        molecule_extents_t extents = molecule_extents_t::molecule_extents_t_gemmi(model, expansion_size);

        EXPECT_NEAR(extents.centre.x(), 26.481f, 1e-4);
        EXPECT_NEAR(extents.centre.y(), 5.563f, 1e-4);
        EXPECT_NEAR(extents.centre.z(), 51.505f, 1e-4);

        EXPECT_NEAR(extents.left.x(), 1.593f, 1e-4);
        EXPECT_NEAR(extents.left.y(), 5.563f, 1e-4);
        EXPECT_NEAR(extents.left.z(), 51.505f, 1e-4);

        EXPECT_NEAR(extents.right.x(), 51.369f, 1e-4);
        EXPECT_NEAR(extents.right.y(), 5.563f, 1e-4);
        EXPECT_NEAR(extents.right.z(), 51.505f, 1e-4);

        EXPECT_NEAR(extents.front.x(), 26.481f, 1e-4);
        EXPECT_NEAR(extents.front.y(), 5.563f, 1e-4);
        EXPECT_NEAR(extents.front.z(), 20.067f, 1e-4);

        EXPECT_NEAR(extents.back.x(), 26.481f, 1e-4);
        EXPECT_NEAR(extents.back.y(), 5.563f, 1e-4);
        EXPECT_NEAR(extents.back.z(), 82.943f, 1e-4);

        EXPECT_NEAR(extents.bottom.x(), 26.481f, 1e-4);
        EXPECT_NEAR(extents.bottom.y(), -14.67f, 1e-4);
        EXPECT_NEAR(extents.bottom.z(), 51.505f, 1e-4);

        EXPECT_NEAR(extents.top.x(), 26.481f, 1e-4);
        EXPECT_NEAR(extents.top.y(), 25.796f, 1e-4);
        EXPECT_NEAR(extents.top.z(), 51.505f, 1e-4);

        EXPECT_NEAR(extents.expansion_size_, 0.0f, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}