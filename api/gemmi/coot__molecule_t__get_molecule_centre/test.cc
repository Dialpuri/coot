#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, GetMoleculeCentre) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    coot::Cartesian c = coot::molecule_t_get_molecule_centre_gemmi(st);
    EXPECT_FLOAT_EQ(c.x(), 29.440775f);
    EXPECT_FLOAT_EQ(c.y(), 5.0443473f);
    EXPECT_FLOAT_EQ(c.z(), 50.240242f);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}