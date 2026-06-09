#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, density_box_t) {
    // Load a PDB to get a valid residue pointer for the parameterized constructor
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Get a residue (first model, first chain, first residue)
    gemmi::Residue *residue_p = nullptr;
    if (!st.models.empty()) {
        gemmi::Model &model = st.models[0];
        if (!model.chains.empty()) {
            gemmi::Chain &chain = model.chains.front();
            if (!chain.residues.empty()) {
                residue_p = &chain.residues.front();
            }
        }
    }

    // Case 1: default constructor — just calls init()
    {
        coot::density_box_t box = coot::density_box_t::density_box_t_gemmi();
        EXPECT_EQ(box.density_box, nullptr);
        EXPECT_EQ(box.residue_p, nullptr);
        EXPECT_EQ(box.n_steps, 0);
        EXPECT_NEAR(box.mean, 0.0, 1e-4);
        EXPECT_NEAR(box.var, -1.0, 1e-4);
        EXPECT_FALSE(box.is_weird);
        EXPECT_TRUE(box.empty());
    }

    // Case 2: parameterized constructor with a small density array and a real residue
    {
        int n_steps = 1;
        int n = 2 * n_steps + 1;  // 3
        int n3 = n * n * n;       // 27
        float *density_box = new float[n3];
        for (int i = 0; i < n3; i++) {
            density_box[i] = static_cast<float>(i * 0.5);
        }

        coot::density_box_t box = coot::density_box_t::density_box_t_gemmi(density_box, residue_p, n_steps);
        EXPECT_NE(box.density_box, nullptr);
        EXPECT_NE(box.residue_p, nullptr);
        EXPECT_EQ(box.n_steps, 1);
        EXPECT_FALSE(box.empty());
        EXPECT_EQ(box.nnn(), 27);
        EXPECT_NEAR(box[0], 0.0f, 1e-4);
        EXPECT_NEAR(box[1], 0.5f, 1e-4);

        delete[] density_box;
    }

    // Case 3: parameterized constructor with null residue (edge case)
    {
        int n_steps = 2;
        int n = 2 * n_steps + 1;  // 5
        int n3 = n * n * n;       // 125
        float *density_box = new float[n3];
        for (int i = 0; i < n3; i++) {
            density_box[i] = 0.0f;
        }

        coot::density_box_t box = coot::density_box_t::density_box_t_gemmi(density_box, nullptr, n_steps);
        EXPECT_NE(box.density_box, nullptr);
        EXPECT_EQ(box.residue_p, nullptr);
        EXPECT_EQ(box.n_steps, 2);
        EXPECT_FALSE(box.empty());
        EXPECT_EQ(box.nnn(), 125);

        delete[] density_box;
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
