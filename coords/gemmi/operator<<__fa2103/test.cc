#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <sstream>
#include "function.hh"

TEST(OracleTest, Cell_Translation) {
    // Case 1: default-constructed (all zeros)
    {
        Cell_Translation ct;
        EXPECT_EQ(ct.us, 0);
        EXPECT_EQ(ct.vs, 0);
        EXPECT_EQ(ct.ws, 0);
        std::ostringstream oss;
        stream_cell_translation_gemmi(oss, ct);
        EXPECT_EQ(oss.str(), "Cell Trans: (0 0 0)");
    }

    // Case 2: constructed with non-zero values
    {
        Cell_Translation ct(1, -2, 3);
        EXPECT_EQ(ct.us, 1);
        EXPECT_EQ(ct.vs, -2);
        EXPECT_EQ(ct.ws, 3);
        std::ostringstream oss;
        stream_cell_translation_gemmi(oss, ct);
        EXPECT_EQ(oss.str(), "Cell Trans: (1 -2 3)");
    }

    // Case 3: negative values
    {
        Cell_Translation ct(-1, -1, -1);
        EXPECT_EQ(ct.us, -1);
        EXPECT_EQ(ct.vs, -1);
        EXPECT_EQ(ct.ws, -1);
        std::ostringstream oss;
        stream_cell_translation_gemmi(oss, ct);
        EXPECT_EQ(oss.str(), "Cell Trans: (-1 -1 -1)");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}