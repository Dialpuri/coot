#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, typed_distances_find_residues_within_window) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Case 1: half_wl = 4
    // In the original, the constructor initializes the map with half_wl=4 (via generate()),
    // then the test calls find_residues_within_window(4) again — appending to existing entries.
    // So we replicate: first call (constructor init), second call (explicit test call).
    {
        std::map<gemmi::CRA, std::vector<gemmi::CRA>, coot::typed_distances::cra_less> residues_within_window;

        // Simulate constructor init: find_residues_within_window(4)
        coot::typed_distances::find_residues_within_window_gemmi(st.models[0], 4, residues_within_window);

        EXPECT_EQ(residues_within_window.size(), 268u);

        // Second explicit call — appends to existing map
        coot::typed_distances::find_residues_within_window_gemmi(st.models[0], 4, residues_within_window);

        EXPECT_EQ(residues_within_window.size(), 268u);

        if (!residues_within_window.empty()) {
            // Find the residue with chain A, seqnum 6, resname LYS
            gemmi::CRA target;
            for (const auto& [cra, neighbors] : residues_within_window) {
                if (cra.chain->name == "A" &&
                    cra.residue->seqid.num.value == 6 &&
                    cra.residue->name == "LYS") {
                    target = cra;
                    break;
                }
            }
            auto it = residues_within_window.find(target);
            EXPECT_TRUE(it != residues_within_window.end());
            if (it != residues_within_window.end()) {
                const gemmi::CRA& res_cra = it->first;

                EXPECT_STREQ(res_cra.chain->name.c_str(), "A");
                EXPECT_EQ(res_cra.residue->seqid.num.value, 6);
                EXPECT_STREQ(res_cra.residue->name.c_str(), "LYS");
                EXPECT_EQ(it->second.size(), 18u);
            }
        }
    }

    // Case 2: edge — half_wl = 0 (only self)
    {
        std::map<gemmi::CRA, std::vector<gemmi::CRA>, coot::typed_distances::cra_less> residues_within_window;
        // In original, td.residues_within_window.clear() happens first
        EXPECT_EQ(residues_within_window.size(), 0u);

        coot::typed_distances::find_residues_within_window_gemmi(st.models[0], 0, residues_within_window);

        EXPECT_EQ(residues_within_window.size(), 268u);

        if (!residues_within_window.empty()) {
            // In half_wl=0, every residue's window has exactly 1 entry (itself)
            auto it2 = residues_within_window.begin();
            EXPECT_EQ(it2->second.size(), 1u);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}