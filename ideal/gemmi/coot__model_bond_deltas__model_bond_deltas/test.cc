#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, model_bond_deltas) {
    // Case 1: default constructor
    {
        coot::model_bond_deltas_gemmi mbd_default;
        EXPECT_TRUE(mbd_default.mol == nullptr);
        EXPECT_TRUE(mbd_default.imol == coot::protein_geometry::IMOL_ENC_UNSET);
        EXPECT_TRUE(mbd_default.geom_p == nullptr);
        EXPECT_EQ(mbd_default.size(), 0u);
    }

    // Case 2: parameterized constructor with real data
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        ASSERT_FALSE(st.models.empty());

        gemmi::Model *mol = &st.models[0];
        ASSERT_NE(mol, nullptr);

        coot::protein_geometry geom;
        geom.init_standard();

        coot::model_bond_deltas_gemmi mbd(mol, 0, &geom);
        EXPECT_NE(mbd.mol, nullptr);
        EXPECT_EQ(mbd.imol, 0);
        EXPECT_NE(mbd.geom_p, nullptr);
        EXPECT_EQ(mbd.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}