#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/math.hpp>
#include <sstream>
#include "function.hh"

TEST(OracleTest, matrix) {
  // Case 1: 3x3 matrix with known values
  {
    gemmi::Mat33 m;
    m.a[0][0] = 1.0; m.a[0][1] = 2.0; m.a[0][2] = 3.0;
    m.a[1][0] = 4.0; m.a[1][1] = 5.0; m.a[1][2] = 6.0;
    m.a[2][0] = 0.0; m.a[2][1] = 0.0; m.a[2][2] = 0.0;

    std::ostringstream oss;
    print_mat33_gemmi(oss, m);
    std::string s = oss.str();
    EXPECT_FALSE(s.empty());
    EXPECT_NE(s.find("1.000000"), std::string::npos);
    EXPECT_NE(s.find("6.000000"), std::string::npos);
  }

  // Case 2: 3x3 matrix with single non-zero value
  {
    gemmi::Mat33 m;
    m.a[0][0] = 42.0;
    m.a[0][1] = 0.0; m.a[0][2] = 0.0;
    m.a[1][0] = 0.0; m.a[1][1] = 0.0; m.a[1][2] = 0.0;
    m.a[2][0] = 0.0; m.a[2][1] = 0.0; m.a[2][2] = 0.0;

    std::ostringstream oss;
    print_mat33_gemmi(oss, m);
    std::string s = oss.str();
    EXPECT_FALSE(s.empty());
    EXPECT_NE(s.find("42.000000"), std::string::npos);
  }

  // Case 3: 3x3 identity matrix
  {
    gemmi::Mat33 m;
    m.a[0][0] = 1.0; m.a[0][1] = 0.0; m.a[0][2] = 0.0;
    m.a[1][0] = 0.0; m.a[1][1] = 1.0; m.a[1][2] = 0.0;
    m.a[2][0] = 0.0; m.a[2][1] = 0.0; m.a[2][2] = 1.0;

    std::ostringstream oss;
    print_mat33_gemmi(oss, m);
    std::string s = oss.str();
    EXPECT_FALSE(s.empty());
    EXPECT_NE(s.find("1.000000"), std::string::npos);
    EXPECT_NE(s.find("0.000000"), std::string::npos);
  }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}