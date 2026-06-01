#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, residues_near_residues_for_residues) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    ASSERT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];

    // Find chain "A"
    gemmi::Chain* chainA = nullptr;
    for (gemmi::Chain& ch : model.chains) {
        if (ch.name == "A") {
            chainA = &ch;
            break;
        }
    }
    ASSERT_NE(chainA, nullptr);

    ASSERT_GE(chainA->residues.size(), 3u);

    // Build all_molecule_map: map some residues to neighbor sets
    std::map<gemmi::Residue*, std::set<gemmi::Residue*>> all_molecule_map;

    std::set<gemmi::Residue*> neighbors;
    neighbors.insert(&chainA->residues[2]);
    neighbors.insert(&chainA->residues[3]);

    all_molecule_map[&chainA->residues[0]] = neighbors;
    all_molecule_map[&chainA->residues[1]] = neighbors;
    all_molecule_map[&chainA->residues[10]] = neighbors;

    // Case 1: limit_to_these_residues_vec contains residues 0 and 1 (but not 10)
    {
        std::set<gemmi::Residue*> limit_set;
        limit_set.insert(&chainA->residues[0]);
        limit_set.insert(&chainA->residues[1]);

        EXPECT_EQ(all_molecule_map.size(), 3u);
        EXPECT_EQ(limit_set.size(), 2u);

        std::map<gemmi::Residue*, std::set<gemmi::Residue*>> result =
            coot::residues_near_residues_for_residues_gemmi(all_molecule_map, limit_set);

        EXPECT_EQ(result.size(), 2u);

        // Print keys in result - iterate and check
        int key_idx = 0;
        for (const auto& kv : result) {
            gemmi::Residue* res = kv.first;
            if (key_idx == 0) {
                EXPECT_EQ(chainA->name, "A");
                EXPECT_EQ(res->seqid.num.value, 1);
                EXPECT_EQ(res->name, "MET");
            } else if (key_idx == 1) {
                EXPECT_EQ(chainA->name, "A");
                EXPECT_EQ(res->seqid.num.value, 2);
                EXPECT_EQ(res->name, "GLU");
            }
            EXPECT_EQ(kv.second.size(), 2u);
            key_idx++;
        }
    }

    // Case 2: empty limit set -> result should be empty
    {
        std::set<gemmi::Residue*> limit_set;
        // leave empty

        EXPECT_EQ(limit_set.size(), 0u);

        std::map<gemmi::Residue*, std::set<gemmi::Residue*>> result =
            coot::residues_near_residues_for_residues_gemmi(all_molecule_map, limit_set);

        EXPECT_EQ(result.size(), 0u);
    }

    // Case 3: limit set contains all keys -> result should equal all_molecule_map size
    {
        std::set<gemmi::Residue*> limit_set;
        limit_set.insert(&chainA->residues[0]);
        limit_set.insert(&chainA->residues[1]);
        limit_set.insert(&chainA->residues[10]);

        EXPECT_EQ(limit_set.size(), 3u);

        std::map<gemmi::Residue*, std::set<gemmi::Residue*>> result =
            coot::residues_near_residues_for_residues_gemmi(all_molecule_map, limit_set);

        EXPECT_EQ(result.size(), 3u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}