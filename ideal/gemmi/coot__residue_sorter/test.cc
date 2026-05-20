#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, ResidueSorter) {
    // Load test data
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty()) << "PDB loading failed";

    // Select all residues in chain A, residues 10-20
    std::vector<gemmi::CRA> selected_cras;
    for (auto& model : st.models)
      for (auto& chain : model.chains)
        if (chain.name == "A")
          for (auto& res : chain.residues) {
            int seq = res.seqid.num.value;
            if (seq >= 10 && seq <= 20) {
              selected_cras.push_back({&chain, &res, nullptr});
            }
          }

    ASSERT_GE(selected_cras.size(), 2u) << "Not enough residues for testing";

    // Test case 1: Direct comparison of adjacent residues (10, 11)
    {
        gemmi::CRA cra10, cra11;
        bool found10 = false, found11 = false;
        for (const auto& c : selected_cras) {
            if (c.residue->seqid.num.value == 10) { cra10 = c; found10 = true; }
            if (c.residue->seqid.num.value == 11) { cra11 = c; found11 = true; }
        }
        ASSERT_TRUE(found10) << "Residue 10 not found";
        ASSERT_TRUE(found11) << "Residue 11 not found";

        std::pair<bool, gemmi::CRA> p10(true, cra10);
        std::pair<bool, gemmi::CRA> p11(true, cra11);

        // residue_sorter(10, 11) should return true (10 comes before 11)
        EXPECT_TRUE(coot::residue_sorter_gemmi(p10, p11));
        // residue_sorter(11, 10) should return false (11 does not come before 10)
        EXPECT_FALSE(coot::residue_sorter_gemmi(p11, p10));
    }

    // Test case 2: Same chain, non-adjacent residues (10, 15)
    {
        gemmi::CRA cra10, cra15;
        bool found10 = false, found15 = false;
        for (const auto& c : selected_cras) {
            if (c.residue->seqid.num.value == 10) { cra10 = c; found10 = true; }
            if (c.residue->seqid.num.value == 15) { cra15 = c; found15 = true; }
        }
        ASSERT_TRUE(found10) << "Residue 10 not found";
        ASSERT_TRUE(found15) << "Residue 15 not found";

        std::pair<bool, gemmi::CRA> p10(true, cra10);
        std::pair<bool, gemmi::CRA> p15(true, cra15);

        // residue_sorter(10, 15) should return true (10 comes before 15)
        EXPECT_TRUE(coot::residue_sorter_gemmi(p10, p15));
        // residue_sorter(15, 10) should return false (15 does not come before 10)
        EXPECT_FALSE(coot::residue_sorter_gemmi(p15, p10));
    }

    // Test case 3: Sorting a vector
    {
        std::vector<std::pair<bool, gemmi::CRA>> residues;
        for (const auto& c : selected_cras)
            residues.push_back(std::make_pair(true, c));

        // Reverse the vector to test sorting
        std::reverse(residues.begin(), residues.end());

        // Check that we have the expected number of residues (11 from 10-20 inclusive)
        EXPECT_EQ(residues.size(), 11u);

        // Sort using residue_sorter_gemmi
        std::sort(residues.begin(), residues.end(), coot::residue_sorter_gemmi);

        // Verify first residue after sorting is residue 10
        EXPECT_EQ(residues[0].second.chain->name, std::string("A"));
        EXPECT_EQ(residues[0].second.residue->seqid.num.value, 10);

        // Verify last residue after sorting is residue 20
        EXPECT_EQ(residues[residues.size()-1].second.chain->name, std::string("A"));
        EXPECT_EQ(residues[residues.size()-1].second.residue->seqid.num.value, 20);

        // Verify sorted - for properly sorted sequence,
        // residue_sorter(prev, next) should be true, or they are equal
        bool sorted = true;
        for (size_t i = 1; i < residues.size(); i++) {
            if (!coot::residue_sorter_gemmi(residues[i-1], residues[i])) {
                // Check if they're equal (both directions return false)
                if (!coot::residue_sorter_gemmi(residues[i], residues[i-1])) {
                    // equal - that's ok
                } else {
                    sorted = false;
                }
            }
            // If residue_sorter(prev, next) is true, that's good (prev < next)
        }
        EXPECT_TRUE(sorted);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}