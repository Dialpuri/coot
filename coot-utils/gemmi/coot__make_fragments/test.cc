#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

static const std::string pdb_file = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";

TEST(OracleTest, make_fragments) {
    // Case 1: Valid PDB with multiple chains and gaps
    {
        gemmi::Structure st = gemmi::read_pdb_file(pdb_file);

        coot::fragment_container_t fc = coot::make_fragments_gemmi(st);

        EXPECT_EQ(fc.ranges.size(), 4);

        // Fragment 0
        {
            const auto &frag = fc.ranges[0];
            EXPECT_EQ(frag.chain_id, "A");
            EXPECT_EQ(frag.residues.size(), 36);
            EXPECT_EQ(frag.start_res.chain_id, "A");
            EXPECT_EQ(frag.start_res.res_no, 1);
            EXPECT_EQ(frag.start_res.ins_code, ' ');
            EXPECT_EQ(frag.end_res.chain_id, "A");
            EXPECT_EQ(frag.end_res.res_no, 36);
            EXPECT_EQ(frag.end_res.ins_code, ' ');
        }

        // Fragment 1
        {
            const auto &frag = fc.ranges[1];
            EXPECT_EQ(frag.chain_id, "A");
            EXPECT_EQ(frag.residues.size(), 26);
            EXPECT_EQ(frag.start_res.chain_id, "A");
            EXPECT_EQ(frag.start_res.res_no, 46);
            EXPECT_EQ(frag.start_res.ins_code, ' ');
            EXPECT_EQ(frag.end_res.chain_id, "A");
            EXPECT_EQ(frag.end_res.res_no, 72);
            EXPECT_EQ(frag.end_res.ins_code, ' ');
        }

        // Fragment 2
        {
            const auto &frag = fc.ranges[2];
            EXPECT_EQ(frag.chain_id, "A");
            EXPECT_EQ(frag.residues.size(), 70);
            EXPECT_EQ(frag.start_res.chain_id, "A");
            EXPECT_EQ(frag.start_res.res_no, 76);
            EXPECT_EQ(frag.start_res.ins_code, ' ');
            EXPECT_EQ(frag.end_res.chain_id, "A");
            EXPECT_EQ(frag.end_res.res_no, 146);
            EXPECT_EQ(frag.end_res.ins_code, ' ');
        }

        // Fragment 3
        {
            const auto &frag = fc.ranges[3];
            EXPECT_EQ(frag.chain_id, "A");
            EXPECT_EQ(frag.residues.size(), 132);
            EXPECT_EQ(frag.start_res.chain_id, "A");
            EXPECT_EQ(frag.start_res.res_no, 166);
            EXPECT_EQ(frag.start_res.ins_code, ' ');
            EXPECT_EQ(frag.end_res.chain_id, "A");
            EXPECT_EQ(frag.end_res.res_no, 298);
            EXPECT_EQ(frag.end_res.ins_code, ' ');
        }
    }

    // Case 2: empty structure — should return empty fragment_container_t
    {
        coot::fragment_container_t fc = coot::make_fragments_gemmi(gemmi::Structure{});
        EXPECT_EQ(fc.ranges.size(), 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}