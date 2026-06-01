#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "mini-mol/mini-mol.hh"
#include "clipper/core/coords.h"

TEST(OracleTest, update_positions_from) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    EXPECT_EQ(st.models.size(), 1u);

    gemmi::Model &model = st.models[0];

    // Get first chain
    ASSERT_GT(model.chains.size(), 0u);
    gemmi::Chain &chain = model.chains[0];

    // Get first residue
    ASSERT_GT(chain.residues.size(), 0u);
    gemmi::Residue &res = chain.residues[0];

    EXPECT_EQ(chain.name, "A");
    EXPECT_EQ(res.seqid.num.value, 1);
    EXPECT_EQ(res.name, "MET");
    EXPECT_EQ(res.atoms.size(), 8u);

    // Helper: build a coot::minimol::residue from a gemmi::Residue
    auto build_minires = [](gemmi::Residue &gres, bool pop_back_one = false) {
        coot::minimol::residue minires;
        for (const gemmi::Atom &gat : gres.atoms) {
            coot::minimol::atom ma;
            float fx = static_cast<float>(gat.pos.x);
            float fy = static_cast<float>(gat.pos.y);
            float fz = static_cast<float>(gat.pos.z);
            ma.pos = clipper::Coord_orth(fx, fy, fz);
            ma.name = gat.name;
            ma.element = gat.element.name();
            minires.atoms.push_back(ma);
        }
        if (pop_back_one) {
            minires.atoms.pop_back();
        }
        minires.seqnum = gres.seqid.num.value;
        return minires;
    };

    // Case 1: same atom count — update_positions_from should work
    {
        coot::minimol::residue minires = build_minires(res);

        // Store original position of first atom
        clipper::Coord_orth original_pos = minires.atoms[0].pos;

        // Overwrite all positions to 999 to prove update restores them
        for (size_t i = 0; i < minires.atoms.size(); i++) {
            minires.atoms[i].pos = clipper::Coord_orth(999.0, 999.0, 999.0);
        }

        clipper::Coord_orth after_overwrite = minires.atoms[0].pos;

        EXPECT_TRUE(minires.atoms.size() == res.atoms.size());
        EXPECT_NEAR(original_pos.x(), 12.334, 1e-4);
        EXPECT_NEAR(original_pos.y(), 2.772, 1e-4);
        EXPECT_NEAR(original_pos.z(), 74.951, 1e-4);
        EXPECT_NEAR(after_overwrite.x(), 999.0, 1e-4);
        EXPECT_NEAR(after_overwrite.y(), 999.0, 1e-4);
        EXPECT_NEAR(after_overwrite.z(), 999.0, 1e-4);

        // Call the function under test
        coot::minimol::update_positions_from_gemmi(minires, &res);

        clipper::Coord_orth after_update = minires.atoms[0].pos;
        EXPECT_NEAR(after_update.x(), 12.334, 1e-4);
        EXPECT_NEAR(after_update.y(), 2.772, 1e-4);
        EXPECT_NEAR(after_update.z(), 74.951, 1e-4);
        // Oracle printed "false" for restored (== comparison)
        EXPECT_FALSE(after_update == original_pos);
    }

    // Case 2: different atom count — guard clause should prevent update
    {
        coot::minimol::residue minires = build_minires(res, true); // pop_back_one

        // Overwrite all positions
        for (size_t i = 0; i < minires.atoms.size(); i++) {
            minires.atoms[i].pos = clipper::Coord_orth(777.0, 777.0, 777.0);
        }

        clipper::Coord_orth before = minires.atoms[0].pos;

        EXPECT_EQ(minires.atoms.size(), 7u);
        EXPECT_EQ(res.atoms.size(), 8u);
        EXPECT_FALSE(minires.atoms.size() == res.atoms.size());

        // Call — should be a no-op due to guard
        coot::minimol::update_positions_from_gemmi(minires, &res);

        clipper::Coord_orth after = minires.atoms[0].pos;
        EXPECT_NEAR(before.x(), 777.0, 1e-4);
        EXPECT_NEAR(before.y(), 777.0, 1e-4);
        EXPECT_NEAR(before.z(), 777.0, 1e-4);
        EXPECT_NEAR(after.x(), 777.0, 1e-4);
        EXPECT_NEAR(after.y(), 777.0, 1e-4);
        EXPECT_NEAR(after.z(), 777.0, 1e-4);
        EXPECT_TRUE(before == after);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}