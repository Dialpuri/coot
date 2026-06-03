#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: find the chain ID for a residue pointer by scanning the model
static std::string chain_id_for_residue(const gemmi::Model& model, const gemmi::Residue* res) {
    for (const auto& chain : model.chains) {
        for (const auto& r : chain.residues) {
            if (&r == res) {
                return chain.name;
            }
        }
    }
    return "";
}

TEST(OracleTest, residues_near_residues) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    ASSERT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];

    // Find residue with seqnum 50 in chain A
    gemmi::CRA target_cra{nullptr, nullptr, nullptr};
    for (auto& chain : model.chains) {
        for (auto& res : chain.residues) {
            if (res.seqid.num.value == 50) {
                target_cra = gemmi::CRA{&chain, &res, nullptr};
                break;
            }
        }
    }
    ASSERT_TRUE(target_cra.residue != nullptr);

    // Case 1: Single residue — find neighbors of //A/50
    {
        std::vector<std::pair<bool, gemmi::CRA>> residues_vec;
        residues_vec.push_back(std::make_pair(true, target_cra));

        float dist_crit = 5.0;

        auto result = coot::residues_near_residues_gemmi(residues_vec, st, dist_crit);

        EXPECT_EQ(result.size(), 1u);

        for (const auto& kv : result) {
            const gemmi::Residue* key = kv.first;
            const std::set<const gemmi::Residue*>& neighbors = kv.second;

            EXPECT_EQ(chain_id_for_residue(model, key), "A");
            EXPECT_EQ(key->seqid.num.value, 50);

            EXPECT_EQ(neighbors.size(), 8u);

            // Verify neighbor identities
            std::set<std::string> neighbor_ids;
            for (const gemmi::Residue* n : neighbors) {
                std::string id = chain_id_for_residue(model, n)
                               + "/" + std::to_string(n->seqid.num.value)
                               + "/" + n->name;
                neighbor_ids.insert(id);
            }
            EXPECT_NE(neighbor_ids.find("A/46/SER"), neighbor_ids.end());
            EXPECT_NE(neighbor_ids.find("A/47/THR"), neighbor_ids.end());
            EXPECT_NE(neighbor_ids.find("A/48/ALA"), neighbor_ids.end());
            EXPECT_NE(neighbor_ids.find("A/49/ILE"), neighbor_ids.end());
            EXPECT_NE(neighbor_ids.find("A/51/GLU"), neighbor_ids.end());
            EXPECT_NE(neighbor_ids.find("A/52/ILE"), neighbor_ids.end());
            EXPECT_NE(neighbor_ids.find("A/53/SER"), neighbor_ids.end());
            EXPECT_NE(neighbor_ids.find("A/54/LEU"), neighbor_ids.end());
        }
    }

    // Case 2: Empty residue vector — should return empty map
    {
        std::vector<std::pair<bool, gemmi::CRA>> residues_vec; // empty
        float dist_crit = 5.0;

        auto result = coot::residues_near_residues_gemmi(residues_vec, st, dist_crit);
        EXPECT_EQ(result.size(), 0u);
    }

    // Case 3: Residue with very small distance — should find fewer neighbors
    {
        std::vector<std::pair<bool, gemmi::CRA>> residues_vec;
        residues_vec.push_back(std::make_pair(true, target_cra));

        float dist_crit = 0.1;

        auto result = coot::residues_near_residues_gemmi(residues_vec, st, dist_crit);
        EXPECT_EQ(result.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}