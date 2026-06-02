#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_lsq_matrix) {
    // Create two gemmi::Structure instances and load the same PDB
    gemmi::Structure st1 = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Structure st2 = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    {
        // Verify model count
        EXPECT_EQ(st1.models.size(), 1);
        EXPECT_EQ(st2.models.size(), 1);
    }

    gemmi::Model* mol1 = &st1.models[0];
    gemmi::Model* mol2 = &st2.models[0];

    // Case 1: Valid match between chain A (ref) and chain B (mov), CA atoms, residues 24-165
    {
        std::vector<coot::lsq_range_match_info_t> matches;
        coot::lsq_range_match_info_t match(
            24,   // to_reference_start_resno
            165,  // to_reference_end_resno
            "A",  // reference_chain_id
            24,   // from_matcher_start_resno
            165,  // from_matcher_end_resno
            "B",  // matcher_chain_id
            2     // match_type_flag: CA
        );
        matches.push_back(match);

        EXPECT_EQ(match.reference_chain_id, "A");
        EXPECT_EQ(match.matcher_chain_id, "B");
        EXPECT_EQ(match.to_reference_start_resno, 24);
        EXPECT_EQ(match.to_reference_end_resno, 165);
        EXPECT_EQ(match.from_matcher_start_resno, 24);
        EXPECT_EQ(match.from_matcher_end_resno, 165);
        EXPECT_EQ(match.match_type_flag, 2);

        int every_nth = 1;
        bool summary_to_screen = false;

        auto result = coot::util::get_lsq_matrix_gemmi(mol1, mol2, matches, every_nth, summary_to_screen);

        EXPECT_EQ(result.first, 0);
        std::string rtop_str = result.second.format();
        EXPECT_FALSE(rtop_str.empty());
        EXPECT_NE(rtop_str.find("0"), std::string::npos);
    }

    // Case 2: Empty matches — should produce istat=0
    {
        std::vector<coot::lsq_range_match_info_t> matches; // empty
        auto result = coot::util::get_lsq_matrix_gemmi(mol1, mol2, matches, 1, false);

        EXPECT_EQ(result.first, 0);
        std::string rtop_str = result.second.format();
        EXPECT_FALSE(rtop_str.empty());
        EXPECT_NE(rtop_str.find("0"), std::string::npos);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}