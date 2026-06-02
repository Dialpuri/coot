#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include <sstream>

TEST(OracleTest, CXXMatrix) {
    // Case 1: identity matrix (default constructor sets diagonal to 1.0)
    {
        matrix4x4 identity;
        std::ostringstream oss;
        print_matrix_gemmi(oss, identity);
        std::string result = oss.str();
        EXPECT_EQ(result, "( [ 1 0 0 0 ]\n [ 0 1 0 0 ]\n [ 0 0 1 0 ]\n [ 0 0 0 1 ])\n");
    }

    // Case 2: matrix from a 4x4 array of doubles
    {
        double m[4][4] = {
            {1, 2, 3, 4},
            {5, 6, 7, 8},
            {9, 10, 11, 12},
            {13, 14, 15, 16}
        };
        matrix4x4 mat(m);
        std::ostringstream oss;
        print_matrix_gemmi(oss, mat);
        std::string result = oss.str();
        EXPECT_EQ(result, "( [ 1 2 3 4 ]\n [ 5 6 7 8 ]\n [ 9 10 11 12 ]\n [ 13 14 15 16 ])\n");
    }

    // Case 3: matrix from quaternion (w=1, x=y=z=0 → identity)
    {
        double quat[] = {0.0, 0.0, 0.0, 1.0}; // x, y, z, w
        matrix4x4 rot = matrix4x4::matrixWithQuaternion(quat);
        std::ostringstream oss;
        print_matrix_gemmi(oss, rot);
        std::string result = oss.str();
        EXPECT_EQ(result, "( [ 1 0 0 0 ]\n [ 0 1 0 0 ]\n [ 0 0 1 0 ]\n [ 0 0 0 1 ])\n");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}