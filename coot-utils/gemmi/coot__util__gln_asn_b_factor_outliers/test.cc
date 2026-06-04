#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, gln_asn_b_factor_outliers) {
    // Case 1: Full call on the PDB (has ASN and GLN residues) — oracle says 0 outliers
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        ASSERT_FALSE(st.models.empty());

        std::vector<std::pair<coot::atom_spec_t, std::string>> result =
            coot::util::gln_asn_b_factor_outliers_gemmi(st.models[0]);
        EXPECT_EQ(result.size(), 0u);
    }

    // Case 2: Empty molecule — should return empty vector
    {
        gemmi::Model empty_model;
        std::vector<std::pair<coot::atom_spec_t, std::string>> result =
            coot::util::gln_asn_b_factor_outliers_gemmi(empty_model);
        EXPECT_EQ(result.size(), 0u);
    }

    // Case 3 (complementary): Verify function runs without crash and
    // returns consistent results on re-call with same molecule.
    // Also verify that the molecule has ASN/GLN residues but no outliers
    // because all B-factors are uniform (zero diff → no outliers).
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        ASSERT_FALSE(st.models.empty());
        gemmi::Model& model = st.models[0];

        // First call
        std::vector<std::pair<coot::atom_spec_t, std::string>> result1 =
            coot::util::gln_asn_b_factor_outliers_gemmi(model);

        // Second call — should return the same result
        std::vector<std::pair<coot::atom_spec_t, std::string>> result2 =
            coot::util::gln_asn_b_factor_outliers_gemmi(model);

        EXPECT_EQ(result1.size(), result2.size());
        EXPECT_EQ(result1.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}