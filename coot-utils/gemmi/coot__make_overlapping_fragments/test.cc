#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, make_overlapping_fragments) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Case 1: valid chain A with fragment_length=19
    {
        std::string chain_id = "A";
        unsigned int fragment_length = 19;
        coot::fragment_container_t fc = coot::make_overlapping_fragments_gemmi(st, chain_id, fragment_length);
        EXPECT_EQ(fc.ranges.size(), 226u);
        if (!fc.ranges.empty()) {
            auto const &r0 = fc.ranges.front();
            EXPECT_EQ(r0.chain_id, "A");
            EXPECT_EQ(r0.start_res.res_no, 1);
            EXPECT_EQ(r0.end_res.res_no, 19);
        }
    }

    // Case 2: non-existent chain — should produce 0 fragments
    {
        std::string chain_id = "Z";
        unsigned int fragment_length = 19;
        coot::fragment_container_t fc = coot::make_overlapping_fragments_gemmi(st, chain_id, fragment_length);
        EXPECT_EQ(fc.ranges.size(), 0u);
    }

    // Case 3: fragment_length=1 — nothing added (≤10 threshold)
    {
        std::string chain_id = "A";
        unsigned int fragment_length = 1;
        coot::fragment_container_t fc = coot::make_overlapping_fragments_gemmi(st, chain_id, fragment_length);
        EXPECT_EQ(fc.ranges.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}