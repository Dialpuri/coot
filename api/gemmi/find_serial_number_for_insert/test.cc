#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, find_serial_number_for_insert) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_EQ(st.models.size(), 1u);
    gemmi::Model &model = st.models[0];

    // Case 1: Insert at seqnum=50 in chain A (middle)
    // Should find first residue with seqnum > 50
    {
        int seqnum_for_new = 50;
        std::string ins_code_for_new = "";
        std::string chain_id = "A";

        std::pair<int, gemmi::Residue *> result =
            find_serial_number_for_insert_gemmi(model, seqnum_for_new, ins_code_for_new, chain_id);

        EXPECT_EQ(result.first, 41);
        EXPECT_NE(result.second, nullptr);
        if (result.second) {
            EXPECT_EQ(result.second->seqid.num.value, 51);
            EXPECT_EQ(std::string(result.second->name), "GLU");
        }
    }

    // Case 2: Insert after all residues (seqnum=999)
    // Should return -1 (no residue with higher seqnum)
    {
        int seqnum_for_new = 999;
        std::string ins_code_for_new = "";
        std::string chain_id = "A";

        std::pair<int, gemmi::Residue *> result =
            find_serial_number_for_insert_gemmi(model, seqnum_for_new, ins_code_for_new, chain_id);

        EXPECT_EQ(result.first, -1);
        EXPECT_EQ(result.second, nullptr);
    }

    // Case 3: Insert before first residue (seqnum=0)
    // Should find residue 1 (serial_no=0)
    {
        int seqnum_for_new = 0;
        std::string ins_code_for_new = "";
        std::string chain_id = "A";

        std::pair<int, gemmi::Residue *> result =
            find_serial_number_for_insert_gemmi(model, seqnum_for_new, ins_code_for_new, chain_id);

        EXPECT_EQ(result.first, 0);
        EXPECT_NE(result.second, nullptr);
        if (result.second) {
            EXPECT_EQ(result.second->seqid.num.value, 1);
            EXPECT_EQ(std::string(result.second->name), "MET");
        }
    }

    // Case 4: Wrong chain (chain X) — should return -1
    {
        int seqnum_for_new = 50;
        std::string ins_code_for_new = "";
        std::string chain_id = "X";

        std::pair<int, gemmi::Residue *> result =
            find_serial_number_for_insert_gemmi(model, seqnum_for_new, ins_code_for_new, chain_id);

        EXPECT_EQ(result.first, -1);
        EXPECT_EQ(result.second, nullptr);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}