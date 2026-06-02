#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include <utility>
#include <sstream>
#include <streambuf>

TEST(OracleTest, optimize_point) {
    // Suppress noise from protein_geometry ("no restraints for type ...")
    std::streambuf *original_cout = std::cout.rdbuf(NULL);

    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    ASSERT_GT(st.models.size(), 0) << "Failed to load PDB";

    gemmi::Model& model = st.models[0];

    clipper::Coord_orth from_pt(10.0, 10.0, 10.0);
    clipper::Coord_orth to_pt(20.0, 20.0, 20.0);

    // Set v_hat manually (normally set by generate())
    clipper::Coord_orth diff = to_pt - from_pt;
    clipper::Coord_orth v_hat(diff.unit());

    // === Case 1: nearby atoms selected using SelectAtoms ===
    {
        clipper::Coord_orth start_pt1(15.0, 15.0, 15.0);
        double radius1 = 7.0;

        // Get atom list for this selection
        std::vector<gemmi::CRA> atom_list1 =
            coot::hole::make_atom_selection_gemmi(model, start_pt1, radius1);

        std::pair<clipper::Coord_orth, double> result1 =
            coot::hole::optimize_point_gemmi(start_pt1, v_hat, model, radius1);

        // Verify oracle output for case 1
        EXPECT_EQ(static_cast<int>(atom_list1.size()), 0);
        EXPECT_EQ(result1.first.x(), 15.0);
        EXPECT_EQ(result1.first.y(), 15.0);
        EXPECT_EQ(result1.first.z(), 15.0);
        EXPECT_EQ(result1.second, -1.0);
    }

    // === Case 2: far from structure (no atoms selected) ===
    {
        clipper::Coord_orth start_pt2(1000.0, 1000.0, 1000.0);
        double radius1 = 7.0;

        // Get atom list for this selection
        std::vector<gemmi::CRA> atom_list2 =
            coot::hole::make_atom_selection_gemmi(model, start_pt2, radius1);

        std::pair<clipper::Coord_orth, double> result2 =
            coot::hole::optimize_point_gemmi(start_pt2, v_hat, model, radius1);

        // Verify oracle output for case 2
        EXPECT_EQ(static_cast<int>(atom_list2.size()), 0);
        EXPECT_EQ(result2.first.x(), 1000.0);
        EXPECT_EQ(result2.first.y(), 1000.0);
        EXPECT_EQ(result2.first.z(), 1000.0);
        EXPECT_EQ(result2.second, -1.0);
    }

    // Restore cout
    std::cout.rdbuf(original_cout);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}