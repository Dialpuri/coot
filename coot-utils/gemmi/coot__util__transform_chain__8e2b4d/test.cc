#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, transform_chain) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_GE(st.models.size(), 1u);

    gemmi::Model &model = st.models[0];
    EXPECT_GE(model.chains.size(), 1);

    gemmi::Chain &chain = model.chains[0];

    // Case 1: Translation by (10, 20, 30) — pure translation
    {
        // Record BEFORE coordinates of first atom
        double bx = 0, by = 0, bz = 0;
        if (!chain.residues.empty() && !chain.residues[0].atoms.empty()) {
            gemmi::Atom &at = chain.residues[0].atoms[0];
            bx = at.pos.x;
            by = at.pos.y;
            bz = at.pos.z;
        }

        EXPECT_NEAR(bx, 12.334, 1e-3);
        EXPECT_NEAR(by, 2.772, 1e-3);
        EXPECT_NEAR(bz, 74.951, 1e-3);

        clipper::Mat33<double> identity(1, 0, 0, 0, 1, 0, 0, 0, 1);
        clipper::Coord_orth translation(10.0, 20.0, 30.0);
        clipper::RTop_orth rtop(identity, translation);

        EXPECT_NO_THROW(coot::util::transform_chain_gemmi(chain, rtop));

        // Record AFTER coordinates of first atom
        if (!chain.residues.empty() && !chain.residues[0].atoms.empty()) {
            gemmi::Atom &at = chain.residues[0].atoms[0];
            EXPECT_NEAR(at.pos.x, 22.334, 1e-4);
            EXPECT_NEAR(at.pos.y, 22.772, 1e-4);
            EXPECT_NEAR(at.pos.z, 104.951, 1e-4);

            // Displacement should be exactly (10, 20, 30)
            EXPECT_NEAR(at.pos.x - bx, 10.0, 1e-4);
            EXPECT_NEAR(at.pos.y - by, 20.0, 1e-4);
            EXPECT_NEAR(at.pos.z - bz, 30.0, 1e-4);
        }
    }

    // Case 2: Identity transform (zero translation) — no-op case
    // NOTE: chain is already translated from Case 1
    {
        double bx = 0, by = 0, bz = 0;
        if (!chain.residues.empty() && !chain.residues[0].atoms.empty()) {
            gemmi::Atom &at = chain.residues[0].atoms[0];
            bx = at.pos.x;
            by = at.pos.y;
            bz = at.pos.z;
        }

        // before = post-case-1 values
        EXPECT_NEAR(bx, 22.334, 1e-3);
        EXPECT_NEAR(by, 22.772, 1e-3);
        EXPECT_NEAR(bz, 104.951, 1e-3);

        clipper::Mat33<double> identity(1, 0, 0, 0, 1, 0, 0, 0, 1);
        clipper::Coord_orth zero_trans(0.0, 0.0, 0.0);
        clipper::RTop_orth rtop(identity, zero_trans);

        EXPECT_NO_THROW(coot::util::transform_chain_gemmi(chain, rtop));

        if (!chain.residues.empty() && !chain.residues[0].atoms.empty()) {
            gemmi::Atom &at = chain.residues[0].atoms[0];
            // Should be unchanged (identity transform)
            EXPECT_NEAR(at.pos.x, 22.334, 1e-4);
            EXPECT_NEAR(at.pos.y, 22.772, 1e-4);
            EXPECT_NEAR(at.pos.z, 104.951, 1e-4);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}