#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, transform_by) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());

    // Get a reference atom to track
    gemmi::Model& model = st.models[0];
    ASSERT_FALSE(model.chains.empty());
    gemmi::Chain& chain = model.chains[0];
    ASSERT_FALSE(chain.residues.empty());
    gemmi::Residue& res = chain.residues[0];
    ASSERT_FALSE(res.atoms.empty());
    gemmi::Atom* atom0 = &res.atoms[0];

    double x_before = atom0->pos.x;
    double y_before = atom0->pos.y;
    double z_before = atom0->pos.z;

    EXPECT_NEAR(x_before, 12.334, 1e-4);
    EXPECT_NEAR(y_before, 2.772, 1e-4);
    EXPECT_NEAR(z_before, 74.951, 1e-4);

    // Case 1: Translate by +10.0 in x
    {
        clipper::Mat33<double> identity(1, 0, 0, 0, 1, 0, 0, 0, 1);
        clipper::Coord_orth translation(10.0, 0.0, 0.0);
        clipper::RTop_orth rtop(identity, translation);

        coot::transform_by_gemmi(st, rtop);

        double x_after = atom0->pos.x;
        double y_after = atom0->pos.y;
        double z_after = atom0->pos.z;

        EXPECT_NEAR(x_after, 22.334, 1e-4);
        EXPECT_NEAR(y_after, 2.772, 1e-4);
        EXPECT_NEAR(z_after, 74.951, 1e-4);
        EXPECT_NEAR(x_after - x_before, 10.0, 1e-4);
        EXPECT_NEAR(y_after - y_before, 0.0, 1e-4);
        EXPECT_NEAR(z_after - z_before, 0.0, 1e-4);
    }

    // Case 2: Identity transform — reload fresh PDB
    gemmi::Structure st2 = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st2.models.empty());

    gemmi::Model& model2 = st2.models[0];
    gemmi::Chain& chain2 = model2.chains[0];
    gemmi::Residue& res2 = chain2.residues[0];
    gemmi::Atom* atom1 = &res2.atoms[0];

    double x_before2 = atom1->pos.x;
    double y_before2 = atom1->pos.y;
    double z_before2 = atom1->pos.z;

    EXPECT_NEAR(x_before2, 12.334, 1e-4);
    EXPECT_NEAR(y_before2, 2.772, 1e-4);
    EXPECT_NEAR(z_before2, 74.951, 1e-4);

    {
        clipper::Mat33<double> identity(1, 0, 0, 0, 1, 0, 0, 0, 1);
        clipper::Coord_orth zero_trans(0.0, 0.0, 0.0);
        clipper::RTop_orth rtop(identity, zero_trans);

        coot::transform_by_gemmi(st2, rtop);

        double x_after2 = atom1->pos.x;
        double y_after2 = atom1->pos.y;
        double z_after2 = atom1->pos.z;

        EXPECT_NEAR(x_after2, 12.334, 1e-4);
        EXPECT_NEAR(y_after2, 2.772, 1e-4);
        EXPECT_NEAR(z_after2, 74.951, 1e-4);
        EXPECT_NEAR(x_after2 - x_before2, 0.0, 1e-4);
        EXPECT_NEAR(y_after2 - y_before2, 0.0, 1e-4);
        EXPECT_NEAR(z_after2 - z_before2, 0.0, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}