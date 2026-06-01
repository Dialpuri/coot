#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_residue) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Case 1: valid residue — chain A, residue 10
    {
        coot::residue_spec_t rs("A", 10);
        gemmi::Residue *res = coot::util::get_residue_gemmi(rs, st);
        EXPECT_NE(res, nullptr);
        EXPECT_EQ(std::string(res->name), "ILE");
        EXPECT_EQ(res->seqid.num.value, 10);
    }

    // Case 2: invalid residue — chain A, residue 9999 (doesn't exist)
    {
        coot::residue_spec_t rs("A", 9999);
        gemmi::Residue *res = coot::util::get_residue_gemmi(rs, st);
        EXPECT_EQ(res, nullptr);
    }

    // Case 3: invalid chain — chain Z (doesn't exist), residue 1
    {
        coot::residue_spec_t rs("Z", 1);
        gemmi::Residue *res = coot::util::get_residue_gemmi(rs, st);
        EXPECT_EQ(res, nullptr);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}