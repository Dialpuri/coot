#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, extents) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Case 1: valid residue spec — chain A, residue 10
    {
        std::vector<coot::residue_spec_t> specs;
        specs.push_back(coot::residue_spec_t("A", 10));

        std::pair<clipper::Coord_orth, clipper::Coord_orth> e =
            coot::util::extents_gemmi(st, specs);

        EXPECT_NEAR(e.first.x(), 26.018, 1e-3);
        EXPECT_NEAR(e.first.y(), 8.313, 1e-3);
        EXPECT_NEAR(e.first.z(), 65.461, 1e-3);
        EXPECT_NEAR(e.second.x(), 28.347, 1e-3);
        EXPECT_NEAR(e.second.y(), 12.742, 1e-3);
        EXPECT_NEAR(e.second.z(), 69.153, 1e-3);
    }

    // Case 2: multiple residue specs — chain A, residues 10, 20, 30
    {
        std::vector<coot::residue_spec_t> specs;
        specs.push_back(coot::residue_spec_t("A", 10));
        specs.push_back(coot::residue_spec_t("A", 20));
        specs.push_back(coot::residue_spec_t("A", 30));

        std::pair<clipper::Coord_orth, clipper::Coord_orth> e =
            coot::util::extents_gemmi(st, specs);

        EXPECT_NEAR(e.first.x(), 22.897, 1e-3);
        EXPECT_NEAR(e.first.y(), -0.246, 1e-3);
        EXPECT_NEAR(e.first.z(), 65.461, 1e-3);
        EXPECT_NEAR(e.second.x(), 28.347, 1e-3);
        EXPECT_NEAR(e.second.y(), 12.742, 1e-3);
        EXPECT_NEAR(e.second.z(), 73.577, 1e-3);
    }

    // Case 3: empty specs vector — edge case, should return default extents
    {
        std::vector<coot::residue_spec_t> specs;

        std::pair<clipper::Coord_orth, clipper::Coord_orth> e =
            coot::util::extents_gemmi(st, specs);

        EXPECT_NEAR(e.first.x(), 99999, 1e-3);
        EXPECT_NEAR(e.first.y(), 99999, 1e-3);
        EXPECT_NEAR(e.first.z(), 99999, 1e-3);
        EXPECT_NEAR(e.second.x(), -99999, 1e-3);
        EXPECT_NEAR(e.second.y(), -99999, 1e-3);
        EXPECT_NEAR(e.second.z(), -99999, 1e-3);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}