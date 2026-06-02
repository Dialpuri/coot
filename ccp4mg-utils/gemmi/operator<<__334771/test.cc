#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <cmath>
#include <sstream>
#include "function.hh"

TEST(OracleTest, Quat) {
    // Case 1: identity quaternion (default constructor)
    {
        Quat q;
        std::ostringstream os;
        print_quat_gemmi(os, q);
        EXPECT_EQ(os.str(), "1 0 0 0");
    }

    // Case 2: axis-angle — 90° about Z axis (wi==1)
    {
        Quat q(0, 0, 1, 1, 90.0);
        std::ostringstream os;
        print_quat_gemmi(os, q);
        EXPECT_EQ(os.str(), "0.707107 0 0 0.707107");
    }

    // Case 3: euler angles (wi==0)
    {
        Quat q(0.0, 0.0, 1.5707963, 0, 0);
        std::ostringstream os;
        print_quat_gemmi(os, q);
        EXPECT_EQ(os.str(), "0.999906 0 0 0.0137074");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}