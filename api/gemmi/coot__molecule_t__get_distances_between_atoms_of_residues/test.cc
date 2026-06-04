#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_distances_between_atoms_of_residues) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty()) << "PDB load failed";

    gemmi::Model &model = st.models[0];

    // Case 1: two adjacent residues with a reasonable distance cutoff
    {
        std::string cid1 = "//A/10";
        std::string cid2 = "//A/11";
        float dist_max = 5.0;
        std::vector<coot::atom_distance_t> results =
            coot::get_distances_between_atoms_of_residues_gemmi(cid1, cid2, dist_max, model);

        EXPECT_EQ(results.size(), 20u);
        if (!results.empty()) {
            EXPECT_EQ(results[0].atom_1.chain_id, "A");
            EXPECT_EQ(results[0].atom_1.res_no, 10);
            EXPECT_EQ(results[0].atom_1.atom_name, " N  ");
            EXPECT_EQ(results[0].atom_2.chain_id, "A");
            EXPECT_EQ(results[0].atom_2.res_no, 11);
            EXPECT_EQ(results[0].atom_2.atom_name, " N  ");
            EXPECT_NEAR(results[0].distance, 2.67663, 1e-4);
        }
    }

    // Case 2: residues on different chains — should return empty
    {
        std::string cid1 = "//A/10";
        std::string cid2 = "//B/10";
        float dist_max = 5.0;
        std::vector<coot::atom_distance_t> results =
            coot::get_distances_between_atoms_of_residues_gemmi(cid1, cid2, dist_max, model);

        EXPECT_EQ(results.size(), 0u);
    }

    // Case 3: non-existent residue — should return empty
    {
        std::string cid1 = "//A/10";
        std::string cid2 = "//A/9999";
        float dist_max = 5.0;
        std::vector<coot::atom_distance_t> results =
            coot::get_distances_between_atoms_of_residues_gemmi(cid1, cid2, dist_max, model);

        EXPECT_EQ(results.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}