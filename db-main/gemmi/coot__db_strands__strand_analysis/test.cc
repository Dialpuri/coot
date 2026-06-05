#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <vector>
#include <string>
#include <cmath>
#include "function.hh"

TEST(strand_analysis_test, strand_analysis) {
    // Load the PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-dev/coot/api/gideondoesntapprove.pdb");

    ASSERT_FALSE(st.models.empty());

    // Check sheets (should already be populated from PDB SHEET records)
    EXPECT_EQ(st.sheets.size(), 4u);

    // We need to pick a reasonable strand_length. Let's inspect one strand
    // from the sheets data to see what lengths exist.
    int max_strand_len = 0;
    for (const auto& sheet : st.sheets) {
        for (const auto& strand : sheet.strands) {
            int len = std::abs(strand.end.res_id.seqid.num.value - strand.start.res_id.seqid.num.value) + 1;
            if (len > max_strand_len) max_strand_len = len;
        }
    }
    EXPECT_EQ(max_strand_len, 7);

    // Case 1: strand_length = max found (should match at least one)
    {
        int strand_length = max_strand_len;
        std::vector<coot::minimol::molecule> rv =
            coot::db_strands::strand_analysis_gemmi(st, strand_length);
        EXPECT_EQ(rv.size(), 2u);
        EXPECT_EQ(rv[0].get_number_of_atoms(), 35);
        EXPECT_EQ(rv[1].get_number_of_atoms(), 35);
    }

    // Case 2: strand_length = 1 (too short, should filter out all)
    {
        int strand_length = 1;
        std::vector<coot::minimol::molecule> rv =
            coot::db_strands::strand_analysis_gemmi(st, strand_length);
        EXPECT_EQ(rv.size(), 0u);
    }

    // Case 3: strand_length very large (no match)
    {
        int strand_length = 9999;
        std::vector<coot::minimol::molecule> rv =
            coot::db_strands::strand_analysis_gemmi(st, strand_length);
        EXPECT_EQ(rv.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}