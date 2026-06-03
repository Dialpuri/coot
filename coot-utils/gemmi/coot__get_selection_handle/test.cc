#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_selection_handle_gemmi) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Case 1: valid atom spec (chain A, residue 10, atom CA)
    {
        coot::atom_spec_t spec;
        spec.chain_id   = "A";
        spec.res_no     = 10;
        spec.ins_code   = "";
        spec.atom_name  = "CA";
        spec.alt_conf   = "";

        int sel = coot::get_selection_handle_gemmi(&st, spec);
        EXPECT_EQ(sel, 2);
        if (sel >= 0) {
            coot::delete_selection_gemmi(sel);
        }
    }

    // Case 2: atom that does not exist (residue 9999)
    {
        coot::atom_spec_t spec;
        spec.chain_id   = "A";
        spec.res_no     = 9999;
        spec.ins_code   = "";
        spec.atom_name  = "CA";
        spec.alt_conf   = "";

        int sel = coot::get_selection_handle_gemmi(&st, spec);
        EXPECT_EQ(sel, 2);
        if (sel >= 0) {
            coot::delete_selection_gemmi(sel);
        }
    }

    // Case 3: nullptr manager — should return -1
    {
        coot::atom_spec_t spec;
        spec.chain_id   = "A";
        spec.res_no     = 10;
        spec.ins_code   = "";
        spec.atom_name  = "CA";
        spec.alt_conf   = "";

        int sel = coot::get_selection_handle_gemmi(nullptr, spec);
        EXPECT_EQ(sel, -1);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}