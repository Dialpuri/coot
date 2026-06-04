#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, dots_representation_info_t) {
    // Case 1: default constructor — dots_representation_info_t()
    {
        pli::dots_representation_info_t dots = pli::dots_representation_info_t::dots_representation_info_t_gemmi();
        EXPECT_TRUE(dots.is_open_p());
        EXPECT_EQ(dots.get_name(), "");
    }

    // Case 2: string constructor — dots_representation_info_t(const std::string &n)
    {
        pli::dots_representation_info_t dots = pli::dots_representation_info_t::dots_representation_info_t_gemmi("my_dots");
        EXPECT_TRUE(dots.is_open_p());
        EXPECT_EQ(dots.get_name(), "my_dots");
    }

    // Case 3: close_yourself() changes is_open from true to false
    {
        pli::dots_representation_info_t dots = pli::dots_representation_info_t::dots_representation_info_t_gemmi();
        EXPECT_TRUE(dots.is_open_p());
        dots.close_yourself();
        EXPECT_FALSE(dots.is_open_p());
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}