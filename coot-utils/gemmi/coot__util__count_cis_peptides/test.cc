#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// ... TEST(...) blocks here ...

TEST(OracleTest, cis_peptides) {
    // Case 1: count cis peptides in full molecule
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        int n_cis = coot::util::count_cis_peptides_gemmi(st);
        EXPECT_EQ(n_cis, 1);
    }

    // Case 2: cis peptide details
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        std::vector<coot::util::cis_peptide_info_t> cis_info =
            coot::util::cis_peptides_info_from_coords_gemmi(st);
        EXPECT_EQ(cis_info.size(), 1u);
        std::string actual = cis_info[0].chain_id_1 + " " + std::to_string(cis_info[0].resno_1)
                           + " - " + cis_info[0].chain_id_2 + " " + std::to_string(cis_info[0].resno_2);
        EXPECT_EQ(actual, "A 253 - A 254");
    }

    // Case 3: empty structure (no atoms) — edge case
    {
        gemmi::Structure empty_st;
        int n_cis = coot::util::count_cis_peptides_gemmi(empty_st);
        EXPECT_EQ(n_cis, 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}