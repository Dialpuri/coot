#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, MoleculeCentre) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty()) << "PDB load failed";

    // case 1: valid molecule
    {
        std::pair<bool, clipper::Coord_orth> result = coot::centre_of_molecule_gemmi(st);
        EXPECT_TRUE(result.first);
        EXPECT_NEAR(result.second.x(), 29.4408, 1e-4);
        EXPECT_NEAR(result.second.y(), 5.04435, 1e-4);
        EXPECT_NEAR(result.second.z(), 50.2402, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}