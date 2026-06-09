#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, residue) {
    // Case 1: default constructor — sets seqnum to mmdb::MinInt4
    {
        coot::minimol::residue r = coot::minimol::residue::residue_gemmi();
        EXPECT_EQ(r.seqnum, -2147483647);
        EXPECT_EQ(r.ins_code, "");
        EXPECT_EQ(r.name, "");
        EXPECT_EQ(r.n_atoms(), 0);
        EXPECT_TRUE(r.is_empty());
    }

    // Case 2: constructor with int seqnum
    {
        coot::minimol::residue r = coot::minimol::residue::residue_gemmi(42);
        EXPECT_EQ(r.seqnum, 42);
        EXPECT_EQ(r.ins_code, "");
        EXPECT_EQ(r.name, "");
        EXPECT_EQ(r.n_atoms(), 0);
    }

    // Case 3: constructor with int + resname
    {
        coot::minimol::residue r = coot::minimol::residue::residue_gemmi(100, "ALA");
        EXPECT_EQ(r.seqnum, 100);
        EXPECT_EQ(r.ins_code, "");
        EXPECT_EQ(r.name, "ALA");
        EXPECT_EQ(r.n_atoms(), 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
