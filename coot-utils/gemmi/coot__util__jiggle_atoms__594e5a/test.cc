#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include "clipper/core/coords.h"
#include <utility>
#include <vector>

TEST(OracleTest, jiggle_atoms) {
  gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
  ASSERT_FALSE(st.models.empty());

  gemmi::Model &model = st.models[0];
  ASSERT_FALSE(model.chains.empty());

  gemmi::Chain &chain = model.chains[0];

  std::vector<gemmi::Atom const *> atoms;
  for (gemmi::Residue &res : chain.residues) {
    for (gemmi::Atom &at : res.atoms) {
      atoms.push_back(&at);
    }
  }

  // Case 1: typical call with annealing_factor=1.0
  {
    float jiggle_trans_scale_factor = 1.0f;
    float annealing_factor = 1.0f;

    // Compute centre point from atoms
    double cx = 0, cy = 0, cz = 0;
    for (auto const *a : atoms) {
      cx += a->pos.x; cy += a->pos.y; cz += a->pos.z;
    }
    cx /= atoms.size(); cy /= atoms.size(); cz /= atoms.size();
    clipper::Coord_orth centre_pt(cx, cy, cz);

    EXPECT_EQ(atoms.size(), 2106u);

    // Centre point is deterministic from the PDB
    EXPECT_NEAR(centre_pt.x(), 29.4402, 1e-3);
    EXPECT_NEAR(centre_pt.y(), 5.04272, 1e-3);
    EXPECT_NEAR(centre_pt.z(), 50.2466, 1e-3);

    auto result = coot::util::jiggle_atoms_gemmi(atoms, centre_pt, jiggle_trans_scale_factor, annealing_factor);

    // Output atom count must match input
    EXPECT_EQ(result.second.size(), 2106u);

    // Original atom0 coords
    gemmi::Atom const *orig = atoms[0];
    EXPECT_NEAR(orig->pos.x, 12.334, 1e-3);
    EXPECT_NEAR(orig->pos.y, 2.772, 1e-3);
    EXPECT_NEAR(orig->pos.z, 74.951, 1e-3);

    // Jiggled atom should differ from original (random, but should have moved)
    // The translation from the oracle was ~2 units, so at least one coord should differ
    double dx = result.second[0].pos.x - orig->pos.x;
    double dy = result.second[0].pos.y - orig->pos.y;
    double dz = result.second[0].pos.z - orig->pos.z;
    EXPECT_GT(dx*dx + dy*dy + dz*dz, 1.0);

    // Translation vector from rtop
    clipper::Coord_orth t(result.first.trn());
    // With jiggle_trans_scale_factor=1.0, translation should be non-trivial
    double t_mag = t.x()*t.x() + t.y()*t.y() + t.z()*t.z();
    EXPECT_GT(t_mag, 1.0);
  }

  // Case 2: annealing_factor=0 (gets clamped to 1.0 internally)
  {
    float jiggle_trans_scale_factor = 0.5f;
    float annealing_factor = 0.0f;

    double cx = 0, cy = 0, cz = 0;
    for (auto const *a : atoms) {
      cx += a->pos.x; cy += a->pos.y; cz += a->pos.z;
    }
    cx /= atoms.size(); cy /= atoms.size(); cz /= atoms.size();
    clipper::Coord_orth centre_pt(cx, cy, cz);

    EXPECT_EQ(atoms.size(), 2106u);

    auto result = coot::util::jiggle_atoms_gemmi(atoms, centre_pt, jiggle_trans_scale_factor, annealing_factor);

    // Output atom count must match input
    EXPECT_EQ(result.second.size(), 2106u);

    // Jiggled atom should differ from original
    gemmi::Atom const *orig = atoms[0];
    double dx = result.second[0].pos.x - orig->pos.x;
    double dy = result.second[0].pos.y - orig->pos.y;
    double dz = result.second[0].pos.z - orig->pos.z;
    EXPECT_GT(dx*dx + dy*dy + dz*dz, 1.0);

    // Translation with scale_factor=0.5
    clipper::Coord_orth t2(result.first.trn());
    double t_mag = t2.x()*t2.x() + t2.y()*t2.y() + t2.z()*t2.z();
    EXPECT_GT(t_mag, 0.01);
  }

  // Case 3: single atom
  {
    std::vector<gemmi::Atom const *> single_atom = {atoms[0]};
    float jiggle_trans_scale_factor = 2.0f;
    float annealing_factor = 0.5f;

    clipper::Coord_orth centre_pt(atoms[0]->pos.x, atoms[0]->pos.y, atoms[0]->pos.z);

    EXPECT_EQ(single_atom.size(), 1u);

    // Centre point is the atom itself
    EXPECT_NEAR(centre_pt.x(), 12.334, 1e-3);
    EXPECT_NEAR(centre_pt.y(), 2.772, 1e-3);
    EXPECT_NEAR(centre_pt.z(), 74.951, 1e-3);

    auto result = coot::util::jiggle_atoms_gemmi(single_atom, centre_pt, jiggle_trans_scale_factor, annealing_factor);

    EXPECT_EQ(result.second.size(), 1u);

    // Original single atom coords
    EXPECT_NEAR(single_atom[0]->pos.x, 12.334, 1e-3);
    EXPECT_NEAR(single_atom[0]->pos.y, 2.772, 1e-3);
    EXPECT_NEAR(single_atom[0]->pos.z, 74.951, 1e-3);

    // With jiggle_trans_scale_factor=2.0, the jiggled atom should be noticeably different
    double dx = result.second[0].pos.x - single_atom[0]->pos.x;
    double dy = result.second[0].pos.y - single_atom[0]->pos.y;
    double dz = result.second[0].pos.z - single_atom[0]->pos.z;
    EXPECT_GT(dx*dx + dy*dy + dz*dz, 1.0);

    clipper::Coord_orth t3(result.first.trn());
    double t_mag = t3.x()*t3.x() + t3.y()*t3.y() + t3.z()*t3.z();
    EXPECT_GT(t_mag, 1.0);
  }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}