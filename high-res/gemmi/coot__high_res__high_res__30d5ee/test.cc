#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/calculate.hpp>
#include "function.hh"

TEST(OracleTest, high_res) {
    // Case 1: iflag = 0
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        int iflag = 0;
        int atom_count = gemmi::count_atom_sites(st);
        EXPECT_EQ(atom_count, 2107);

        coot::high_res::high_res_gemmi hr(st, iflag);
        EXPECT_EQ(gemmi::count_atom_sites(hr.globular_molecule), 2107);
    }

    // Case 2: iflag = 1
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        int iflag = 1;
        int atom_count = gemmi::count_atom_sites(st);
        EXPECT_EQ(atom_count, 2107);

        coot::high_res::high_res_gemmi hr(st, iflag);
        EXPECT_EQ(gemmi::count_atom_sites(hr.globular_molecule), 2107);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}