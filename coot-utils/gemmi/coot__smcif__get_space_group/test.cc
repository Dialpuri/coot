#include <gtest/gtest.h>
#include <clipper/core/spacegroup.h>
#include "function.hh"

TEST(OracleTest, GetSpaceGroup) {
    // Case 1: Valid single symmetry operation (P1 spacegroup)
    {
        std::vector<std::string> symm_strings = {"x,y,z"};
        std::pair<bool, clipper::Spacegroup> result = coot::smcif::get_space_group_gemmi(symm_strings);
        EXPECT_EQ(result.first, true);
        EXPECT_EQ(result.second.symbol_xhm(), "P 1");
    }
    
    // Case 2: Valid multiple symmetry operations (P-1 spacegroup)
    {
        std::vector<std::string> symm_strings = {"x,y,z", "-x,-y,-z"};
        std::pair<bool, clipper::Spacegroup> result = coot::smcif::get_space_group_gemmi(symm_strings);
        EXPECT_EQ(result.first, true);
        EXPECT_EQ(result.second.symbol_xhm(), "P -1");
    }
    
    // Case 3: Valid but more complete P1 with all 3 symmetry operations
    {
        std::vector<std::string> symm_strings = {"x,y,z", "-x,-y,-z", "-x,y,-z"};
        std::pair<bool, clipper::Spacegroup> result = coot::smcif::get_space_group_gemmi(symm_strings);
        EXPECT_EQ(result.first, true);
        EXPECT_EQ(result.second.symbol_xhm(), "P 1 2/m 1");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}