#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// ... TEST(...) blocks here ...

TEST(OracleTest, dots_representation_info_t) {
    // Case 1: construct with a named string
    {
        std::string label = "TestDots";
        pli::dots_representation_info_t dots = pli::dots_representation_info_t::dots_representation_info_t_gemmi(label);
        EXPECT_EQ(dots.get_name(), "TestDots");
        EXPECT_TRUE(dots.is_open_p());
    }

    // Case 2: construct with empty string
    {
        std::string label = "";
        pli::dots_representation_info_t dots = pli::dots_representation_info_t::dots_representation_info_t_gemmi(label);
        EXPECT_EQ(dots.get_name(), "");
        EXPECT_TRUE(dots.is_open_p());
    }

    // Case 3: construct with special characters
    {
        std::string label = "LZA_1299_A";
        pli::dots_representation_info_t dots = pli::dots_representation_info_t::dots_representation_info_t_gemmi(label);
        EXPECT_EQ(dots.get_name(), "LZA_1299_A");
        EXPECT_TRUE(dots.is_open_p());
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}