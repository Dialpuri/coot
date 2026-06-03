#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, residues_near_residues) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_GE(st.models.size(), 1u);

    // Case 1: reasonable distance criterion (should find neighbors)
    {
        float dist_crit = 5.0f;

        auto result = coot::residues_near_residues_gemmi(&st, dist_crit);

        EXPECT_EQ(result.size(), 268u);

        // Build a lookup by (chain, seqnum, resname) for easier checking
        // We need resname which isn't in residue_spec_t, so iterate the model
        std::map<std::string, const std::set<coot::residue_spec_t>*> lookup;
        for (auto const &pair : result) {
            auto const& key_spec = pair.first;
            // Find matching residue to get resname
            for (auto const &chain : st.models[0].chains) {
                if (chain.name != key_spec.chain_id) continue;
                for (auto const &res : chain.residues) {
                    if (res.seqid.num.value == key_spec.res_no) {
                        std::string key_str = key_spec.chain_id + "," + std::to_string(key_spec.res_no) + "," + res.name;
                        lookup[key_str] = &(pair.second);
                        break;
                    }
                }
            }
        }

        // Check A 223 ASP: 9 neighbors
        auto it = lookup.find("A,223,ASP");
        if (it == lookup.end()) {
            FAIL() << "A,223,ASP not found in map";
        } else {
            EXPECT_EQ(it->second->size(), 9u);
        }

        // Check A 224 GLU: 10 neighbors
        it = lookup.find("A,224,GLU");
        if (it == lookup.end()) {
            FAIL() << "A,224,GLU not found in map";
        } else {
            EXPECT_EQ(it->second->size(), 10u);
        }

        // Check A 225 VAL: 5 neighbors
        it = lookup.find("A,225,VAL");
        if (it == lookup.end()) {
            FAIL() << "A,225,VAL not found in map";
        } else {
            EXPECT_EQ(it->second->size(), 5u);
        }
    }

    // Case 2: very small distance criterion (should find nothing)
    {
        float dist_crit = 0.001f;

        auto result = coot::residues_near_residues_gemmi(&st, dist_crit);

        EXPECT_EQ(result.size(), 0u);
    }

    // Case 3: null mol (should return empty map)
    {
        float dist_crit = 5.0f;

        auto result = coot::residues_near_residues_gemmi(nullptr, dist_crit);

        EXPECT_EQ(result.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}