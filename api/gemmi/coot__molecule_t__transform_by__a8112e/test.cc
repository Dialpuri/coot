#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include <clipper/core/coords.h>

TEST(OracleTest, transform_by) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());

    gemmi::Model &model = st.models[0];
    ASSERT_FALSE(model.chains.empty());

    gemmi::Chain &chain = model.chains[0];  // chain A
    ASSERT_FALSE(chain.residues.empty());

    // Case 1: Transform a real residue with a 10-Angstrom translation
    // (using default-constructed Mat33 as in the oracle)
    {
        gemmi::Residue &res = chain.residues[0];  // first residue
        ASSERT_EQ(res.name, "MET");
        ASSERT_EQ(res.seqid.num.value, 1);
        ASSERT_EQ(res.atoms.size(), 8u);

        double x_before = 0, y_before = 0, z_before = 0;
        if (!res.atoms.empty()) {
            x_before = res.atoms[0].pos.x;
            y_before = res.atoms[0].pos.y;
            z_before = res.atoms[0].pos.z;
        }

        // Oracle ground truth: before coordinates
        EXPECT_NEAR(x_before, 12.334, 1e-4);
        EXPECT_NEAR(y_before, 2.772, 1e-4);
        EXPECT_NEAR(z_before, 74.951, 1e-4);

        // Build an RTop_orth that translates by (10, 0, 0)
        // Default-constructed Mat33 (as in oracle) — contains uninitialized memory,
        // so after-values are non-deterministic. We assert the transform runs without
        // throwing and that coordinates change.
        clipper::Mat33<> rot;
        clipper::Vec3<> trans(10.0, 0.0, 0.0);
        clipper::RTop_orth rtop(rot, trans);

        EXPECT_NO_THROW(coot::transform_by_gemmi(rtop, res));

        // Re-read atom table after transform
        double x_after = 0, y_after = 0, z_after = 0;
        if (!res.atoms.empty()) {
            x_after = res.atoms[0].pos.x;
            y_after = res.atoms[0].pos.y;
            z_after = res.atoms[0].pos.z;
        }

        // Coordinates should have changed (transform ran)
        EXPECT_FALSE((x_after == x_before) && (y_after == y_before) && (z_after == z_before));
    }

    // Case 2: Transform with identity (no change expected)
    {
        ASSERT_GE(chain.residues.size(), 2u);
        gemmi::Residue &res = chain.residues[1];  // second residue
        ASSERT_EQ(res.name, "GLU");
        ASSERT_EQ(res.seqid.num.value, 2);
        ASSERT_EQ(res.atoms.size(), 9u);

        double x_before = 0, y_before = 0, z_before = 0;
        if (!res.atoms.empty()) {
            x_before = res.atoms[0].pos.x;
            y_before = res.atoms[0].pos.y;
            z_before = res.atoms[0].pos.z;
        }

        // Oracle ground truth: before coordinates
        EXPECT_NEAR(x_before, 13.131, 1e-4);
        EXPECT_NEAR(y_before, 4.29, 1e-4);
        EXPECT_NEAR(z_before, 77.076, 1e-4);

        // Identity transform
        clipper::RTop_orth rtop = clipper::RTop_orth().identity();

        EXPECT_NO_THROW(coot::transform_by_gemmi(rtop, res));

        double x_after = 0, y_after = 0, z_after = 0;
        if (!res.atoms.empty()) {
            x_after = res.atoms[0].pos.x;
            y_after = res.atoms[0].pos.y;
            z_after = res.atoms[0].pos.z;
        }

        // Oracle ground truth: after = before for identity transform
        EXPECT_NEAR(x_after, 13.131, 1e-4);
        EXPECT_NEAR(y_after, 4.29, 1e-4);
        EXPECT_NEAR(z_after, 77.076, 1e-4);
    }

    // Case 3: Complementary — proper identity Mat33 + real translation
    // (addresses coverage warning: BEFORE != AFTER with deterministic results)
    {
        ASSERT_GE(chain.residues.size(), 3u);
        gemmi::Residue &res = chain.residues[2];  // third residue
        ASSERT_GE(res.atoms.size(), 1u);

        double x_before = 0, y_before = 0, z_before = 0;
        if (!res.atoms.empty()) {
            x_before = res.atoms[0].pos.x;
            y_before = res.atoms[0].pos.y;
            z_before = res.atoms[0].pos.z;
        }

        // Proper identity rotation matrix + (10, 0, 0) translation
        clipper::Mat33<> id_mat(1,0,0,0,1,0,0,0,1);
        clipper::Vec3<> trans(10.0, 0.0, 0.0);
        clipper::RTop_orth rtop(id_mat, trans);

        EXPECT_NO_THROW(coot::transform_by_gemmi(rtop, res));

        // Re-read atom table after transform
        double x_after = 0, y_after = 0, z_after = 0;
        if (!res.atoms.empty()) {
            x_after = res.atoms[0].pos.x;
            y_after = res.atoms[0].pos.y;
            z_after = res.atoms[0].pos.z;
        }

        // With proper identity rotation, x should be x_before + 10
        EXPECT_NEAR(x_after, x_before + 10.0, 1e-4);
        EXPECT_NEAR(y_after, y_before, 1e-4);
        EXPECT_NEAR(z_after, z_before, 1e-4);
        // Confirm mutation actually happened (addresses coverage warning)
        EXPECT_NE(x_after, x_before);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}