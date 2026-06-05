#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <sstream>
#include "function.hh"

TEST(OracleTest, SymmMatrix) {
    // case: zeroed
    {
        SymmMatrix m;
        for (int i=0;i<4;i++) for (int j=0;j<4;j++) m.mat[i][j] = 0.0;
        std::stringstream s;
        operator_stream_insert_gemmi(s, m);
        std::string result = s.str();
        EXPECT_EQ(result, "0 0 0 0 \n0 0 0 0 \n0 0 0 0 \n0 0 0 0 \n");
    }

    // case: identity
    {
        SymmMatrix m;
        for (int i=0;i<4;i++) for (int j=0;j<4;j++) m.mat[i][j] = (i==j) ? 1.0 : 0.0;
        std::stringstream s;
        operator_stream_insert_gemmi(s, m);
        std::string result = s.str();
        EXPECT_EQ(result, "1 0 0 0 \n0 1 0 0 \n0 0 1 0 \n0 0 0 1 \n");
    }

    // case: sequential_1_to_16
    {
        SymmMatrix m;
        double vals[4][4] = {
            {1.0, 2.0, 3.0, 4.0},
            {5.0, 6.0, 7.0, 8.0},
            {9.0, 10.0, 11.0, 12.0},
            {13.0, 14.0, 15.0, 16.0}
        };
        for (int i=0;i<4;i++) for (int j=0;j<4;j++) m.mat[i][j] = vals[i][j];
        std::stringstream s;
        operator_stream_insert_gemmi(s, m);
        std::string result = s.str();
        EXPECT_EQ(result, "1 2 3 4 \n5 6 7 8 \n9 10 11 12 \n13 14 15 16 \n");
    }

    // case: negative_and_fractional (complements the oracle's coverage)
    {
        SymmMatrix m;
        double vals[4][4] = {
            {-1.5, 0.5, -2.25, 0.0},
            {3.33, -0.01, 4.44, -5.55},
            {0.0, 0.0, 0.0, 0.0},
            {-1.0, -2.0, -3.0, -4.0}
        };
        for (int i=0;i<4;i++) for (int j=0;j<4;j++) m.mat[i][j] = vals[i][j];
        std::stringstream s;
        operator_stream_insert_gemmi(s, m);
        std::string result = s.str();
        EXPECT_FALSE(result.empty());
        EXPECT_NE(result.find("-1.5"), std::string::npos);
        EXPECT_NE(result.find("0.5"), std::string::npos);
        EXPECT_NE(result.find("-2.25"), std::string::npos);
        EXPECT_NE(result.find("3.33"), std::string::npos);
        EXPECT_NE(result.find("-4"), std::string::npos);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}