#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include <sstream>
#include <iostream>

TEST(OracleTest, molecule_extents_t) {
  gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
  ASSERT_FALSE(st.models.empty());
  gemmi::Model& model = st.models[0];

  // Case 1: default padding (0)
  {
    coot::molecule_extents_t_gemmi extents = coot::compute_extents_gemmi(model, 0.0);
    std::ostringstream oss;
    oss << std::setprecision(6);  // default stream precision
    coot::operator_stream_insert_gemmi(oss, extents);
    std::string result = oss.str();

    EXPECT_NE(result.find("front:  (26.481,5.563,20.067)"), std::string::npos);
    EXPECT_NE(result.find("back :  (26.481,5.563,82.943)"), std::string::npos);
    EXPECT_NE(result.find("left :  (1.593,5.563,51.505)"), std::string::npos);
    EXPECT_NE(result.find("right:  (51.369,5.563,51.505)"), std::string::npos);
    EXPECT_NE(result.find("top  :  (26.481,25.796,51.505)"), std::string::npos);
    EXPECT_NE(result.find("bottom: (26.481,-14.67,51.505)"), std::string::npos);
  }

  // Case 2: non-zero padding
  {
    coot::molecule_extents_t_gemmi extents = coot::compute_extents_gemmi(model, 5.0);
    std::ostringstream oss;
    oss << std::setprecision(6);
    coot::operator_stream_insert_gemmi(oss, extents);
    std::string result = oss.str();

    EXPECT_NE(result.find("front:  (26.481,5.563,15.067)"), std::string::npos);
    EXPECT_NE(result.find("back :  (26.481,5.563,87.943)"), std::string::npos);
    EXPECT_NE(result.find("left :  (-3.407,5.563,51.505)"), std::string::npos);
    EXPECT_NE(result.find("right:  (56.369,5.563,51.505)"), std::string::npos);
    EXPECT_NE(result.find("top  :  (26.481,30.796,51.505)"), std::string::npos);
    EXPECT_NE(result.find("bottom: (26.481,-19.67,51.505)"), std::string::npos);
  }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}