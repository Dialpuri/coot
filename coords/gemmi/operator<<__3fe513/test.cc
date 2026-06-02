#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/math.hpp>
#include <sstream>
#include "function.hh"

TEST(OracleTest, coot_mat44) {
    // Case 1: identity matrix
    {
        gemmi::Transform t;
        t.set_identity();

        EXPECT_FALSE(coot::is_close_to_unit_matrix_transform_gemmi(t));

        std::ostringstream oss;
        coot::operator_stream_insert_gemmi(oss, t);
        std::string output = oss.str();

        EXPECT_NE(output.find("coot::mat_44"), std::string::npos);
        EXPECT_NE(output.find("[ 1 0 0 0 ]"), std::string::npos);
        EXPECT_NE(output.find("[ 0 1 0 0 ]"), std::string::npos);
        EXPECT_NE(output.find("[ 0 0 1 0 ]"), std::string::npos);
        EXPECT_NE(output.find("[ 0 0 0 1 ]"), std::string::npos);
    }

    // Case 2: translation matrix (translate by 5.0 in x)
    {
        gemmi::Transform t;
        t.set_identity();
        t.vec.x = 5.0;

        EXPECT_FALSE(coot::is_close_to_unit_matrix_transform_gemmi(t));

        std::ostringstream oss;
        coot::operator_stream_insert_gemmi(oss, t);
        std::string output = oss.str();

        EXPECT_NE(output.find("coot::mat_44"), std::string::npos);
        EXPECT_NE(output.find("[ 1 0 0 5 ]"), std::string::npos);
        EXPECT_NE(output.find("[ 0 1 0 0 ]"), std::string::npos);
        EXPECT_NE(output.find("[ 0 0 1 0 ]"), std::string::npos);
        EXPECT_NE(output.find("[ 0 0 0 1 ]"), std::string::npos);
    }

    // Case 3: zero matrix
    {
        gemmi::Transform t;
        // Explicitly zero out the matrix and translation
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                t.mat.a[i][j] = 0.0;
        t.vec.x = 0.0;
        t.vec.y = 0.0;
        t.vec.z = 0.0;

        EXPECT_FALSE(coot::is_close_to_unit_matrix_transform_gemmi(t));

        std::ostringstream oss;
        coot::operator_stream_insert_gemmi(oss, t);
        std::string output = oss.str();

        EXPECT_NE(output.find("coot::mat_44"), std::string::npos);
        EXPECT_NE(output.find("[ 0 0 0 0 ]"), std::string::npos);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}