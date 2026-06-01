#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, residue_spec_t_get_residue) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Case 1: valid residue (chain A, residue 50 = ARG, has atoms)
    {
        coot::residue_spec_t spec("A", 50);
        EXPECT_EQ(spec.chain_id, "A");
        EXPECT_EQ(spec.res_no, 50);
        gemmi::Residue *res = coot::get_residue_gemmi(spec, &st);
        EXPECT_NE(res, nullptr);
        if (res) {
            EXPECT_EQ(std::string(res->name), "ARG");
            EXPECT_EQ(res->seqid.num.value, 50);
            EXPECT_EQ(res->atoms.size(), 11);
        }
    }

    // Case 2: nonexistent residue number
    {
        coot::residue_spec_t spec("A", 9999);
        EXPECT_EQ(spec.chain_id, "A");
        EXPECT_EQ(spec.res_no, 9999);
        gemmi::Residue *res = coot::get_residue_gemmi(spec, &st);
        EXPECT_EQ(res, nullptr);
    }

    // Case 3: nonexistent chain
    {
        coot::residue_spec_t spec("Z", 50);
        EXPECT_EQ(spec.chain_id, "Z");
        EXPECT_EQ(spec.res_no, 50);
        gemmi::Residue *res = coot::get_residue_gemmi(spec, &st);
        EXPECT_EQ(res, nullptr);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}