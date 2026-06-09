#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, is_member_p) {
    std::vector<std::string> v = {"ALA", "GLY", "VAL", "LEU"};

    // Case 1: element present in vector
    {
        std::string a = "GLY";
        bool result = coot::is_member_p_gemmi(v, a);
        EXPECT_EQ(result, true);
    }

    // Case 2: element not present in vector
    {
        std::string a = "SER";
        bool result = coot::is_member_p_gemmi(v, a);
        EXPECT_EQ(result, false);
    }

    // Case 3: empty vector
    {
        std::vector<std::string> empty;
        std::string a = "ALA";
        bool result = coot::is_member_p_gemmi(empty, a);
        EXPECT_EQ(result, false);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
