#include <gtest/gtest.h>
#include "function.hh"

TEST(OracleTest, make_rtop_orth_from) {
    // Test case 1: Identity matrix (no transformation)
    {
        double mat[4][4];
        // Initialize to identity
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                mat[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }

        clipper::RTop_orth rtop = coot::util::make_rtop_orth_from_gemmi(mat);

        clipper::Mat33<double> rot = rtop.rot();
        EXPECT_NEAR(rot(0,0), 1.0, 1e-4);
        EXPECT_NEAR(rot(0,1), 0.0, 1e-4);
        EXPECT_NEAR(rot(0,2), 0.0, 1e-4);
        EXPECT_NEAR(rot(1,0), 0.0, 1e-4);
        EXPECT_NEAR(rot(1,1), 1.0, 1e-4);
        EXPECT_NEAR(rot(1,2), 0.0, 1e-4);
        EXPECT_NEAR(rot(2,0), 0.0, 1e-4);
        EXPECT_NEAR(rot(2,1), 0.0, 1e-4);
        EXPECT_NEAR(rot(2,2), 1.0, 1e-4);

        EXPECT_NEAR(rtop.trn()[0], 0.0, 1e-4);
        EXPECT_NEAR(rtop.trn()[1], 0.0, 1e-4);
        EXPECT_NEAR(rtop.trn()[2], 0.0, 1e-4);
    }

    // Test case 2: Simple translation
    {
        double mat[4][4];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                mat[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }
        // Add translation
        mat[0][3] = 10.5;
        mat[1][3] = -5.25;
        mat[2][3] = 3.0;

        clipper::RTop_orth rtop = coot::util::make_rtop_orth_from_gemmi(mat);

        clipper::Mat33<double> rot = rtop.rot();
        EXPECT_NEAR(rot(0,0), 1.0, 1e-4);
        EXPECT_NEAR(rot(0,1), 0.0, 1e-4);
        EXPECT_NEAR(rot(0,2), 0.0, 1e-4);
        EXPECT_NEAR(rot(1,0), 0.0, 1e-4);
        EXPECT_NEAR(rot(1,1), 1.0, 1e-4);
        EXPECT_NEAR(rot(1,2), 0.0, 1e-4);
        EXPECT_NEAR(rot(2,0), 0.0, 1e-4);
        EXPECT_NEAR(rot(2,1), 0.0, 1e-4);
        EXPECT_NEAR(rot(2,2), 1.0, 1e-4);

        EXPECT_NEAR(rtop.trn()[0], 10.5, 1e-4);
        EXPECT_NEAR(rtop.trn()[1], -5.25, 1e-4);
        EXPECT_NEAR(rtop.trn()[2], 3.0, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}