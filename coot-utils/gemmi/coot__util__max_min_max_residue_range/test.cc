#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Case 1: valid structure
TEST(OracleTest, max_min_max_residue_range_valid_structure) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    int result = coot::util::max_min_max_residue_range_gemmi(st);
    EXPECT_EQ(result, 298);
}

// Case 2: empty structure (no models)
TEST(OracleTest, max_min_max_residue_range_empty_structure) {
    gemmi::Structure st;
    int result = coot::util::max_min_max_residue_range_gemmi(st);
    EXPECT_EQ(result, -1);
}