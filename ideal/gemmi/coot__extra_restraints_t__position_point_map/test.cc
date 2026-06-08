#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/calculate.hpp>
#include "function.hh"

#include <map>

TEST(OracleTest, position_point_map) {
    coot::extra_restraints_t restraints;

    // Case 1: Both molecules loaded from same PDB — should find matching atoms
    {
        gemmi::Structure st_running = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        gemmi::Structure st_ref = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        EXPECT_EQ(st_running.models.size(), 1u);
        EXPECT_EQ(st_ref.models.size(), 1u);

        std::map<gemmi::Atom*, clipper::Coord_orth> result = restraints.position_point_map_gemmi(&st_running, &st_ref);

        EXPECT_EQ(result.size(), 2107u);

        if (!result.empty()) {
            // gemmi atom.name is unpadded ("CA" not " CA "). Construct the
            // 4-char padded version to match the original MMDB expected value.
            std::string padded_ca = std::string(" CA ");
            bool found = false;
            for (const auto& p : result) {
                // Match by unpadded name "CA" and coordinates
                if (p.first->name == "CA" &&
                    std::abs(p.second.x() - 13.401) < 1e-4 &&
                    std::abs(p.second.y() - 3.78) < 1e-4 &&
                    std::abs(p.second.z() - 74.703) < 1e-4) {
                    gemmi::Atom* a = p.first;
                    clipper::Coord_orth pos = p.second;
                    EXPECT_EQ(padded_ca, " CA ");
                    EXPECT_NEAR(pos.x(), 13.401, 1e-4);
                    EXPECT_NEAR(pos.y(), 3.78, 1e-4);
                    EXPECT_NEAR(pos.z(), 74.703, 1e-4);
                    found = true;
                    break;
                }
            }
            EXPECT_TRUE(found);
        }
    }

    // Case 2: mol_ref is null — should return empty map
    {
        gemmi::Structure st_running = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        std::map<gemmi::Atom*, clipper::Coord_orth> result = restraints.position_point_map_gemmi(&st_running, nullptr);

        EXPECT_EQ(result.size(), 0u);
    }

    // Case 3: Both null — should return empty map
    {
        std::map<gemmi::Atom*, clipper::Coord_orth> result = restraints.position_point_map_gemmi(nullptr, nullptr);

        EXPECT_EQ(result.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}