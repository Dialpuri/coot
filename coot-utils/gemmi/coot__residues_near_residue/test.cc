#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include "geometry/residue-and-atom-specs.hh"
#include "function.hh"

TEST(OracleTest, residuesNearResidue) {
    // Load PDB file
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty()) << "Failed to load PDB";

    // Case 1: Valid residue with nearby residues
    {
        std::string cid = "//A/10";
        coot::residue_spec_t rs = coot::cid_to_residue_spec(cid);
        std::vector<coot::residue_spec_t> result = coot::residues_near_residue_gemmi(rs, st, 5.0);
        EXPECT_EQ(result.size(), 9);
    }

    // Case 2: Invalid CID (nonexistent residue)
    {
        std::string cid = "//A/9999";
        coot::residue_spec_t rs = coot::cid_to_residue_spec(cid);
        std::vector<coot::residue_spec_t> result = coot::residues_near_residue_gemmi(rs, st, 5.0);
        EXPECT_EQ(result.size(), 0);
    }

    // Case 3: Valid residue with smaller radius
    {
        std::string cid = "//A/20";
        coot::residue_spec_t rs = coot::cid_to_residue_spec(cid);
        std::vector<coot::residue_spec_t> result = coot::residues_near_residue_gemmi(rs, st, 3.0);
        EXPECT_EQ(result.size(), 4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}