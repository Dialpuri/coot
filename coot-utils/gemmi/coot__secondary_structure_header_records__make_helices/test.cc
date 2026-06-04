#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, make_helices) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    ASSERT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];

    // Case 1: Valid helix with real residues
    {
        int helices_before = static_cast<int>(st.helices.size());
        EXPECT_EQ(helices_before, 0);

        ASSERT_FALSE(model.chains.empty());
        gemmi::Chain& chain = model.chains[0];

        ASSERT_FALSE(chain.residues.empty());
        gemmi::Residue& start_res = chain.residues[0];
        ASSERT_LT(4, chain.residues.size());
        gemmi::Residue& end_res = chain.residues[4];

        // Verify oracle input values
        EXPECT_EQ(chain.name, "A");
        EXPECT_EQ(start_res.seqid.num.value, 1);
        EXPECT_EQ(start_res.name, "MET");
        EXPECT_EQ(chain.name, "A");
        EXPECT_EQ(end_res.seqid.num.value, 5);
        EXPECT_EQ(end_res.name, "GLN");

        coot::secondary_structure_header_records::helix_info_t h_info(
            gemmi::CRA{&chain, &start_res, nullptr},
            gemmi::CRA{&chain, &end_res, nullptr},
            5);
        std::vector<coot::secondary_structure_header_records::helix_info_t> helices;
        helices.push_back(h_info);

        // Verify helix_info_t properties (intermediate observation for coverage)
        EXPECT_EQ(h_info.length, 5);

        coot::secondary_structure_header_records::make_helices_gemmi(st, helices);

        int helices_after = static_cast<int>(st.helices.size());
        EXPECT_EQ(helices_after, 1);

        // Additional coverage: inspect the created helix
        const gemmi::Helix& helix = st.helices[0];
        EXPECT_EQ(helix.start.chain_name, "A");
        EXPECT_EQ(helix.length, 5);
    }

    // Case 2: Empty helix list — should be a no-op
    {
        int helices_before_case2 = static_cast<int>(st.helices.size());
        EXPECT_EQ(helices_before_case2, 1);

        std::vector<coot::secondary_structure_header_records::helix_info_t> empty_helices;
        coot::secondary_structure_header_records::make_helices_gemmi(st, empty_helices);

        int helices_after_case2 = static_cast<int>(st.helices.size());
        EXPECT_EQ(helices_after_case2, 1);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}