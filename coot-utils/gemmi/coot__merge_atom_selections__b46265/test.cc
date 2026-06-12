#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/to_pdb.hpp>
#include "function.hh"

#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

TEST(OracleTest, merge_atom_selections) {

    // --- Case 1: Duplicate chain with overlap → should merge ---
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        int n_models_before = static_cast<int>(st.models.size());
        EXPECT_EQ(n_models_before, 1);

        gemmi::Model* model = &st.models[0];
        int n_chains_before = 0;
        if (model) {
            n_chains_before = static_cast<int>(model->chains.size());
        }
        EXPECT_EQ(n_chains_before, 2);

        std::string chains_before;
        if (model) {
            for (int i = 0; i < n_chains_before; i++) {
                if (!chains_before.empty()) chains_before += ", ";
                chains_before += model->chains[i].name;
            }
        }
        EXPECT_EQ(chains_before, "A, B");

        // Create modified PDB with overlapping chains using gemmi APIs
        // Remove original chain B (HOH), then duplicate chain A as chain B with incremented residues
        {
            gemmi::Structure dup_st;
            dup_st.cell = st.cell;
            dup_st.spacegroup_hm = st.spacegroup_hm;
            dup_st.models.emplace_back();

            // Copy only chain A
            for (gemmi::Chain& chain : st.models[0].chains) {
                if (chain.name == "A") {
                    dup_st.models[0].chains.push_back(chain);
                }
            }

            // Now duplicate chain A as chain B with incremented residues
            gemmi::Chain& chain_a = dup_st.models[0].chains[0];
            gemmi::Chain chain_b = chain_a;
            chain_b.name = "B";
            for (gemmi::Residue& res : chain_b.residues) {
                res.seqid.num.value += 1;
            }
            dup_st.models[0].chains.push_back(chain_b);

            std::ofstream outfile("/tmp/oracle-merge-test.pdb");
            gemmi::write_pdb(dup_st, outfile);
            outfile.close();
        }

        st = gemmi::read_pdb_file("/tmp/oracle-merge-test.pdb");

        model = &st.models[0];
        int n_chains_input = 0;
        if (model) {
            n_chains_input = static_cast<int>(model->chains.size());
        }
        EXPECT_EQ(n_chains_input, 2);

        std::string chains_before_merge;
        if (model) {
            for (int i = 0; i < n_chains_input; i++) {
                if (!chains_before_merge.empty()) chains_before_merge += ", ";
                chains_before_merge += model->chains[i].name;
            }
        }
        EXPECT_EQ(chains_before_merge, "A, B");

        // Call merge
        coot::merge_atom_selections_gemmi(st);

        model = &st.models[0];
        int n_chains_after = 0;
        if (model) {
            n_chains_after = static_cast<int>(model->chains.size());
        }
        EXPECT_EQ(n_chains_after, 1);

        std::string chains_after_merge;
        if (model) {
            for (int i = 0; i < n_chains_after; i++) {
                if (!chains_after_merge.empty()) chains_after_merge += ", ";
                chains_after_merge += model->chains[i].name;
            }
        }
        EXPECT_EQ(chains_after_merge, "B");
        EXPECT_TRUE(n_chains_after < n_chains_input);
    }

    // --- Case 2: No overlapping chains (original PDB) → should NOT merge ---
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        gemmi::Model* model = &st.models[0];
        int n_chains_before = 0;
        if (model) {
            n_chains_before = static_cast<int>(model->chains.size());
        }
        EXPECT_EQ(n_chains_before, 2);

        // Call merge on original PDB — chains A (protein) and B (HOH) should NOT merge
        coot::merge_atom_selections_gemmi(st);

        model = &st.models[0];
        int n_chains_after = 0;
        if (model) {
            n_chains_after = static_cast<int>(model->chains.size());
        }
        EXPECT_EQ(n_chains_after, 2);
        EXPECT_FALSE(n_chains_after < n_chains_before);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}