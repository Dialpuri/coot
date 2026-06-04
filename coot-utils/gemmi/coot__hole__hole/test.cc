#include <gtest/gtest.h>
#include <gemmi/mmread.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(HoleGemmiTest, hole_construction) {
    // Load the CIF
    gemmi::Structure st = gemmi::read_structure_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");

    // Init protein geometry (REQUIRED)
    coot::protein_geometry geom;
    geom.init_standard();

    // Get the gemmi::Model& from the structure
    gemmi::Model& model = st.models[0];

    // --- Case 1: Typical construction with real points ---
    {
        clipper::Coord_orth from_pt(10.0, 20.0, 30.0);
        clipper::Coord_orth to_pt(50.0, 60.0, 70.0);

        auto h = coot::hole::hole_gemmi(model, from_pt, to_pt, geom);

        EXPECT_EQ(h.model, &model);
        EXPECT_NEAR(h.from_pt.x(), 10.0, 1e-4);
        EXPECT_NEAR(h.from_pt.y(), 20.0, 1e-4);
        EXPECT_NEAR(h.from_pt.z(), 30.0, 1e-4);
        EXPECT_NEAR(h.to_pt.x(), 50.0, 1e-4);
        EXPECT_NEAR(h.to_pt.y(), 60.0, 1e-4);
        EXPECT_NEAR(h.to_pt.z(), 70.0, 1e-4);
        EXPECT_EQ(h.colour_map_multiplier, 1);
        EXPECT_EQ(h.colour_map_offset, 0);
        EXPECT_EQ(h.radius_handle, 16777217);
    }

    // --- Case 2: Edge case — from == to (zero-length vector) ---
    {
        clipper::Coord_orth from_pt(0.0, 0.0, 0.0);
        clipper::Coord_orth to_pt(0.0, 0.0, 0.0);

        auto h = coot::hole::hole_gemmi(model, from_pt, to_pt, geom);

        EXPECT_EQ(h.model, &model);
        EXPECT_NEAR(h.from_pt.x(), 0.0, 1e-4);
        EXPECT_NEAR(h.from_pt.y(), 0.0, 1e-4);
        EXPECT_NEAR(h.from_pt.z(), 0.0, 1e-4);
        EXPECT_NEAR(h.to_pt.x(), 0.0, 1e-4);
        EXPECT_NEAR(h.to_pt.y(), 0.0, 1e-4);
        EXPECT_NEAR(h.to_pt.z(), 0.0, 1e-4);
        EXPECT_EQ(h.colour_map_multiplier, 1);
        EXPECT_EQ(h.colour_map_offset, 0);
        EXPECT_EQ(h.radius_handle, 16777217);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}