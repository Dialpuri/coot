#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, shift_to_origin) {
    // Case 1: PDB with manually set crystallography
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        // Manually set cell and spacegroup since the PDB has no CRYST1
        st.cell = gemmi::UnitCell(50.0, 60.0, 70.0, 90.0, 90.0, 90.0);
        st.spacegroup_hm = "P 1";

        clipper::Coord_frac result = coot::util::shift_to_origin_gemmi(st);
        EXPECT_NEAR(result.u(), -1, 1e-4);
        EXPECT_NEAR(result.v(), -0, 1e-4);
        EXPECT_NEAR(result.w(), -1, 1e-4);
    }

    // Case 2: PDB without crystallography
    {
        gemmi::Structure st2 = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        // Do NOT set cell or spacegroup

        clipper::Coord_frac result2 = coot::util::shift_to_origin_gemmi(st2);
        EXPECT_NEAR(result2.u(), -1, 1e-4);
        EXPECT_NEAR(result2.v(), -0, 1e-4);
        EXPECT_NEAR(result2.w(), -1, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}