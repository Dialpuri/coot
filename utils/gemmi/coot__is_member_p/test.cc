#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, is_member_p) {
    // Test case 1: element is in the vector
    {
        std::vector<std::string> v = {"ALA", "GLY", "SER"};
        std::string a = "GLY";
        bool result = coot::is_member_p_gemmi(v, a);
        EXPECT_TRUE(result);
    }

    // Test case 2: element is NOT in the vector
    {
        std::vector<std::string> v = {"ALA", "GLY", "SER"};
        std::string a = "ASP";
        bool result = coot::is_member_p_gemmi(v, a);
        EXPECT_FALSE(result);
    }

    // Test case 3: empty vector
    {
        std::vector<std::string> v;
        std::string a = "ALA";
        bool result = coot::is_member_p_gemmi(v, a);
        EXPECT_FALSE(result);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}