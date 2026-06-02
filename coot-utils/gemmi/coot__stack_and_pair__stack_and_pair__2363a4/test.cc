#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: find a residue by chain ID and sequence number in a gemmi Model
static const gemmi::Residue* find_residue(
    const gemmi::Model& model, char chain_id, int seq_num) {
    for (const auto& chain : model.chains) {
        if (chain.name == std::string(1, chain_id)) {
            for (const auto& res : chain.residues) {
                if (res.seqid.num.value == seq_num) {
                    return &res;
                }
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, stack_and_pair) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    const gemmi::Model& model = st.models[0];

    // Case 1: empty residues_vec — tests guard clause in constructor
    {
        std::vector<std::pair<bool, const gemmi::Residue*>> residues_vec;
        auto result = coot::stack_and_pair::stack_and_pair_gemmi(residues_vec);

        EXPECT_EQ(residues_vec.size(), 0u);
        EXPECT_EQ(result.normal_map.size(), 0u);
    }

    // Case 2: one protein residue — tests non-empty but non-nucleic acid
    {
        const gemmi::Residue* res = find_residue(model, 'A', 10);
        std::vector<std::pair<bool, const gemmi::Residue*>> residues_vec;
        if (res) {
            residues_vec.push_back(std::make_pair(false, res));
        }
        auto result = coot::stack_and_pair::stack_and_pair_gemmi(residues_vec);

        std::string res_info = res ? std::string(res->name) : "nullptr";
        EXPECT_EQ(residues_vec.size(), 1u);
        EXPECT_EQ(res_info, "ILE");
        EXPECT_EQ(result.normal_map.size(), 0u);
    }

    // Case 3: several protein residues — tests multiple entries
    {
        std::vector<std::pair<bool, const gemmi::Residue*>> residues_vec;
        for (int seq : {10, 11, 12}) {
            const gemmi::Residue* res = find_residue(model, 'A', seq);
            if (res) {
                residues_vec.push_back(std::make_pair(false, res));
            }
        }
        auto result = coot::stack_and_pair::stack_and_pair_gemmi(residues_vec);

        EXPECT_EQ(residues_vec.size(), 3u);
        EXPECT_EQ(result.normal_map.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}