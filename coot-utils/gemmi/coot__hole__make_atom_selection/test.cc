#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTestGemmi, hole_make_atom_selection) {
  gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

  ASSERT_EQ(st.models.size(), 1u);
  gemmi::Model& model = st.models[0];

  // Case 1: Center of protein, reasonable radius
  {
    clipper::Coord_orth pt(15.0, 15.0, 15.0);
    double radius_prev = 20.0;

    auto sel = coot::hole::make_atom_selection_gemmi(model, pt, radius_prev);

    EXPECT_EQ(sel.size(), 157u);
  }

  // Case 2: Far from protein, small radius — expect 0 atoms
  {
    clipper::Coord_orth pt(1000.0, 1000.0, 1000.0);
    double radius_prev = 5.0;

    auto sel = coot::hole::make_atom_selection_gemmi(model, pt, radius_prev);

    EXPECT_EQ(sel.size(), 0u);
  }

  // Case 3: Very large radius covering entire protein
  {
    clipper::Coord_orth pt(15.0, 15.0, 15.0);
    double radius_prev = 100.0;

    auto sel = coot::hole::make_atom_selection_gemmi(model, pt, radius_prev);

    EXPECT_EQ(sel.size(), 2106u);
  }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}