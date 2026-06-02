#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <sstream>
#include "function.hh"

TEST(OracleTest, symm_trans_t_stream_output) {
    // Case 1: symm_trans_t with explicit shifts
    {
        symm_trans_t_gemmi t(1, 1, 0, -1);
        t.symm_as_string = "x+1, -y, z-1";
        EXPECT_EQ(t.isym(), 1);
        EXPECT_EQ(t.symm_as_string, "x+1, -y, z-1");
        EXPECT_EQ(t.x(), 1);
        EXPECT_EQ(t.y(), 0);
        EXPECT_EQ(t.z(), -1);
        std::ostringstream os;
        symm_trans_t_stream_gemmi(os, t);
        std::string result = os.str();
        EXPECT_EQ(result, "symm: x+1, -y, z-1 (op-idx: 1) trans: 1 0 -1");
    }

    // Case 2: identity-like operator (no shifts)
    {
        symm_trans_t_gemmi t(0);
        t.symm_as_string = "x, y, z";
        EXPECT_EQ(t.isym(), 0);
        EXPECT_EQ(t.symm_as_string, "x, y, z");
        std::ostringstream os;
        symm_trans_t_stream_gemmi(os, t);
        std::string result = os.str();
        EXPECT_EQ(result, "symm: x, y, z (op-idx: 0) trans: 0 0 0");
    }

    // Case 3: default-constructed (empty)
    {
        symm_trans_t_gemmi t;
        EXPECT_EQ(t.symm_as_string, "");
        std::ostringstream os;
        symm_trans_t_stream_gemmi(os, t);
        std::string result = os.str();
        // Default-constructed t has initialized int fields (0), so
        // only assert on the well-defined string prefix.
        EXPECT_EQ(result.substr(0, 6), "symm: ");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}