#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_residue) {
    // Load the PDB into a gemmi::Structure
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Case 1: Valid residue — chain A, residue 10
    {
        coot::residue_spec_t spec;
        spec.chain_id = "A";
        spec.res_no   = 10;
        spec.ins_code = "";

        gemmi::Residue *res = coot::lidia_utils::get_residue_gemmi(spec, &st);

        EXPECT_NE(res, nullptr);
        if (res) {
            EXPECT_EQ(res->name, "ILE");
            EXPECT_EQ(res->seqid.num.value, 10);
        }
    }

    // Case 2: Invalid residue number — no such residue
    {
        coot::residue_spec_t spec;
        spec.chain_id = "A";
        spec.res_no   = 9999;
        spec.ins_code = "";

        gemmi::Residue *res = coot::lidia_utils::get_residue_gemmi(spec, &st);

        EXPECT_EQ(res, nullptr);
    }

    // Case 3: Nonexistent chain
    {
        coot::residue_spec_t spec;
        spec.chain_id = "Z";
        spec.res_no   = 1;
        spec.ins_code = "";

        gemmi::Residue *res = coot::lidia_utils::get_residue_gemmi(spec, &st);

        EXPECT_EQ(res, nullptr);
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}