#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, chain_atoms_segid) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model& model = st.models[0];

    int n_chains = static_cast<int>(model.chains.size());
    EXPECT_EQ(n_chains, 2);

    // Helper: find a chain by ID
    auto find_chain = [&](const std::string& id) -> const gemmi::Chain* {
        for (auto& c : model.chains) {
            if (c.name == id) return &c;
        }
        return nullptr;
    };

    // Case 1: Chain A — empty segid
    {
        const gemmi::Chain* chain_p = find_chain("A");
        ASSERT_NE(chain_p, nullptr);
        std::string segid = coot::chain_atoms_segid_gemmi(chain_p);
        EXPECT_EQ(segid, "");
    }

    // Case 2: Chain B — empty segid
    {
        const gemmi::Chain* chain_p = find_chain("B");
        ASSERT_NE(chain_p, nullptr);
        std::string segid = coot::chain_atoms_segid_gemmi(chain_p);
        EXPECT_EQ(segid, "");
    }

    // Case 3: No empty chain exists in example.pdb
    {
        bool found_empty = false;
        for (auto& c : model.chains) {
            if (c.residues.empty()) {
                found_empty = true;
                break;
            }
        }
        EXPECT_FALSE(found_empty);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}