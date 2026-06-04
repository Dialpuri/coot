#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/math.hpp>
#include "function.hh"

TEST(OracleTest, matrix_convert) {
    // Case 1: identity rotation + translation (5,10,15)
    {
        gemmi::Mat33 rot(
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0);
        gemmi::Vec3 trans(5.0, 10.0, 15.0);
        gemmi::Transform transform(rot, trans);

        clipper::RTop_orth rtop = coot::util::matrix_convert_gemmi(transform);

        EXPECT_DOUBLE_EQ(rot[0][0], 1.0);
        EXPECT_DOUBLE_EQ(trans.x, 5.0);
        EXPECT_DOUBLE_EQ(trans.y, 10.0);
        EXPECT_DOUBLE_EQ(trans.z, 15.0);
        EXPECT_EQ(rtop.trn().format(), "(         5,        10,        15)");
        EXPECT_EQ(
            rtop.rot().format(),
            "|         1,         0,         0|\n"
            "|         0,         1,         0|\n"
            "|         0,         0,         1|");
    }

    // Case 2: 90° rotation about Z + zero translation
    {
        gemmi::Mat33 rot(
            0.0,  -1.0, 0.0,
            1.0,   0.0, 0.0,
            0.0,   0.0, 1.0);
        gemmi::Vec3 trans(0.0, 0.0, 0.0);
        gemmi::Transform transform(rot, trans);

        clipper::RTop_orth rtop = coot::util::matrix_convert_gemmi(transform);

        EXPECT_DOUBLE_EQ(rot[0][0], 0.0);
        EXPECT_DOUBLE_EQ(rot[0][1], -1.0);
        EXPECT_DOUBLE_EQ(rot[1][0], 1.0);
        EXPECT_DOUBLE_EQ(rot[1][1], 0.0);
        EXPECT_DOUBLE_EQ(trans.x, 0.0);
        EXPECT_EQ(rtop.trn().format(), "(         0,         0,         0)");
        EXPECT_EQ(
            rtop.rot().format(),
            "|         0,        -1,         0|\n"
            "|         1,         0,         0|\n"
            "|         0,         0,         1|");
    }

    // Case 3: all zeros (edge case — degenerate matrix)
    {
        gemmi::Mat33 rot(
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0);
        gemmi::Vec3 trans(0.0, 0.0, 0.0);
        gemmi::Transform transform(rot, trans);

        clipper::RTop_orth rtop = coot::util::matrix_convert_gemmi(transform);

        EXPECT_DOUBLE_EQ(rot[0][0], 0.0);
        EXPECT_EQ(rtop.trn().format(), "(         0,         0,         0)");
        EXPECT_EQ(
            rtop.rot().format(),
            "|         0,         0,         0|\n"
            "|         0,         0,         0|\n"
            "|         0,         0,         0|");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}