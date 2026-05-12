#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <clipper/clipper.h>
#include "function.hh"

TEST(OracleTest, get_residue_closest_to) {
    {
        // case: first test (point near residue 237)
        coot::molecule_t mol;
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        clipper::Coord_orth pt(10.0, 20.0, 30.0);
        coot::residue_spec_t spec1 = mol.get_residue_closest_to_gemmi(st, pt);
        EXPECT_EQ(spec1.model_number, 1);
        EXPECT_EQ(spec1.chain_id, "A");
        EXPECT_EQ(spec1.res_no, 237);
        EXPECT_EQ(spec1.ins_code, "");
    }
    {
        // case: second test (point outside, closest residue 96)
        coot::molecule_t mol;
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        clipper::Coord_orth pt2(100.0, 100.0, 100.0);
        coot::residue_spec_t spec2 = mol.get_residue_closest_to_gemmi(st, pt2);
        EXPECT_EQ(spec2.model_number, 1);
        EXPECT_EQ(spec2.chain_id, "A");
        EXPECT_EQ(spec2.res_no, 96);
        EXPECT_EQ(spec2.ins_code, "");
    }
    {
        // case: third test (different molecule manager, residue 237)
        coot::molecule_t mol;
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        clipper::Coord_orth pt3(15.0, 25.0, 35.0);
        coot::residue_spec_t spec3 = mol.get_residue_closest_to_gemmi(st, pt3);
        EXPECT_EQ(spec3.model_number, 1);
        EXPECT_EQ(spec3.chain_id, "A");
        EXPECT_EQ(spec3.res_no, 237);
        EXPECT_EQ(spec3.ins_code, "");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}