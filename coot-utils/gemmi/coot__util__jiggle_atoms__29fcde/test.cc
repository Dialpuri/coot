#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, jiggle_atoms) {
    // Load the PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    ASSERT_FALSE(st.models.empty());

    gemmi::Model* model = &st.models[0];
    ASSERT_FALSE(model->chains.empty());

    // Gather atoms from first chain
    gemmi::Chain* chain = &model->chains[0];

    std::vector<gemmi::Atom> atoms;
    double sum_x = 0, sum_y = 0, sum_z = 0;

    for (gemmi::Residue& res : chain->residues) {
        for (gemmi::Atom& a : res.atoms) {
            gemmi::Atom copy = a;
            atoms.push_back(copy);
            sum_x += a.pos.x;
            sum_y += a.pos.y;
            sum_z += a.pos.z;
        }
    }

    ASSERT_FALSE(atoms.empty());

    double inv_n = 1.0 / double(atoms.size());
    clipper::Coord_orth centre_pt(sum_x * inv_n, sum_y * inv_n, sum_z * inv_n);

    // ---- Case 1: normal jiggle with scale_factor=1.0, annealing=1.0 ----
    {
        float jiggle_scale_factor = 1.0f;
        float annealing_factor = 1.0f;

        EXPECT_EQ(atoms.size(), 2106u);
        EXPECT_NEAR(centre_pt.x(), 29.4402, 1e-4);
        EXPECT_NEAR(centre_pt.y(), 5.04272, 1e-4);
        EXPECT_NEAR(centre_pt.z(), 50.2466, 1e-4);

        auto result = coot::util::jiggle_atoms_gemmi(atoms, centre_pt, jiggle_scale_factor, annealing_factor);

        EXPECT_EQ(result.second.size(), 2106u);
        EXPECT_NEAR(atoms[0].pos.x, 12.334, 1e-4);
        EXPECT_NEAR(atoms[0].pos.y, 2.772, 1e-4);
        EXPECT_NEAR(atoms[0].pos.z, 74.951, 1e-4);

        // Jiggled position must differ from original (random rotation applied)
        EXPECT_FALSE(
            (std::abs(result.second[0].pos.x - atoms[0].pos.x) < 1e-6) &&
            (std::abs(result.second[0].pos.y - atoms[0].pos.y) < 1e-6) &&
            (std::abs(result.second[0].pos.z - atoms[0].pos.z) < 1e-6)
        ) << "jiggled position should differ from original";
    }

    // ---- Case 2: zero jiggle scale factor (should produce zero translation) ----
    {
        float jiggle_scale_factor = 0.0f;
        float annealing_factor = 1.0f;

        EXPECT_EQ(atoms.size(), 2106u);

        auto result = coot::util::jiggle_atoms_gemmi(atoms, centre_pt, jiggle_scale_factor, annealing_factor);

        EXPECT_EQ(result.second.size(), 2106u);
        EXPECT_NEAR(atoms[0].pos.x, 12.334, 1e-4);
        EXPECT_NEAR(atoms[0].pos.y, 2.772, 1e-4);
        EXPECT_NEAR(atoms[0].pos.z, 74.951, 1e-4);

        // With scale=0, translation should be (0,0,0)
        clipper::Coord_orth t2(result.first.trn());
        EXPECT_NEAR(t2.x(), 0, 1e-4);
        EXPECT_NEAR(t2.y(), 0, 1e-4);
        EXPECT_NEAR(t2.z(), 0, 1e-4);
    }

    // ---- Case 3: large scale factor with annealing=0 (max jiggle; annealing clamped to 1.0) ----
    {
        float jiggle_scale_factor = 5.0f;
        float annealing_factor = 0.0f;

        EXPECT_EQ(atoms.size(), 2106u);

        auto result = coot::util::jiggle_atoms_gemmi(atoms, centre_pt, jiggle_scale_factor, annealing_factor);

        EXPECT_EQ(result.second.size(), 2106u);
        EXPECT_NEAR(atoms[0].pos.x, 12.334, 1e-4);
        EXPECT_NEAR(atoms[0].pos.y, 2.772, 1e-4);
        EXPECT_NEAR(atoms[0].pos.z, 74.951, 1e-4);

        // Jiggled position must differ from original (random rotation applied)
        EXPECT_FALSE(
            (std::abs(result.second[0].pos.x - atoms[0].pos.x) < 1e-6) &&
            (std::abs(result.second[0].pos.y - atoms[0].pos.y) < 1e-6) &&
            (std::abs(result.second[0].pos.z - atoms[0].pos.z) < 1e-6)
        ) << "jiggled position should differ from original";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}