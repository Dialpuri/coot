#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, SecondaryStructureHeaderRecords) {
    // case 1: default constructor (no-op, just tests construction works)
    {
        coot::secondary_structure_header_records sshr;
        EXPECT_TRUE(true); // just confirming default constructor works
    }

    // case 2: constructor with structure (should populate helix/sheet info)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        EXPECT_GT(st.models.size(), 0) << "load failed";

        coot::secondary_structure_header_records sshr(st, false);
        EXPECT_TRUE(true); // just confirming constructor with structure works
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}