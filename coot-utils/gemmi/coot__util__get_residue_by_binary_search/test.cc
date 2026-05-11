#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, GetResidueByBinarySearch) {
    // Setup: Load structure exactly as in oracle
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    { // case 1: existing residue
        std::string chain_id = "A";
        int res_no = 10;
        std::string ins_code = " ";
        gemmi::Residue *res = coot::util::get_residue_by_binary_search_gemmi(chain_id, res_no, ins_code, st);
        EXPECT_TRUE(res != nullptr) << "Residue should be found";
        if (res) {
            EXPECT_EQ(res->seqid.num.value, res_no);
            EXPECT_EQ(std::string(1, res->seqid.icode), ins_code);
            EXPECT_EQ(res->name, "ILE");
        }
    }

    { // case 2: non-existent residue
        std::string chain_id = "A";
        int res_no = 9999;
        std::string ins_code = "";
        gemmi::Residue *res = coot::util::get_residue_by_binary_search_gemmi(chain_id, res_no, ins_code, st);
        EXPECT_TRUE(res == nullptr) << "Residue should not be found";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}