#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, atom_indices_in_other_molecule) {
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Helper to collect CRA for a given chain/residue pair.
    auto select_cras = [&](const std::string &chain_id, int resno) {
        std::vector<gemmi::CRA> result;
        for (auto &model : st.models) {
            for (auto &chain : model.chains) {
                if (chain.name != chain_id) continue;
                for (auto &residue : chain.residues) {
                    if (residue.seqid.num.value != resno) continue;
                    for (auto &atom : residue.atoms) {
                        result.push_back({&chain, &residue, &atom});
                    }
                }
            }
        }
        return result;
    };

    // --- Case 1: identical selections (chain A, res 1) ---
    {
        auto sel1 = select_cras("A", 1);
        auto sel2 = select_cras("A", 1);

        std::set<int> result =
            coot::atom_indices_in_other_molecule_gemmi(sel1, sel2);

        EXPECT_EQ(static_cast<int>(sel1.size()), 8);
        EXPECT_EQ(static_cast<int>(sel2.size()), 8);
        EXPECT_EQ(static_cast<int>(result.size()), 8);
        int i = 0;
        for (auto idx : result) {
            if (i < 10) {
                EXPECT_EQ(idx, i);
                i++;
            }
        }
    }

    // --- Case 2: different residues (no overlap expected) ---
    {
        auto sel1 = select_cras("A", 1);
        auto sel2 = select_cras("A", 2);

        std::set<int> result =
            coot::atom_indices_in_other_molecule_gemmi(sel1, sel2);

        EXPECT_EQ(static_cast<int>(sel1.size()), 8);
        EXPECT_EQ(static_cast<int>(sel2.size()), 9);
        EXPECT_EQ(static_cast<int>(result.size()), 0);
    }

    // --- Case 3: different chains (no overlap expected) ---
    {
        auto sel1 = select_cras("A", 1);
        auto sel2 = select_cras("B", 1);

        std::set<int> result =
            coot::atom_indices_in_other_molecule_gemmi(sel1, sel2);

        EXPECT_EQ(static_cast<int>(sel1.size()), 8);
        EXPECT_EQ(static_cast<int>(sel2.size()), 1);
        EXPECT_EQ(static_cast<int>(result.size()), 0);
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}