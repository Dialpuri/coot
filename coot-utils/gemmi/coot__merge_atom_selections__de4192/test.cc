#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include <fstream>
#include <cstdio>

// Helper: build a selection vector from a gemmi::Model for a given chain name
static std::vector<std::pair<gemmi::Atom*, gemmi::Residue*>>
build_selection(const gemmi::Model& model, const std::string& chain_name) {
    std::vector<std::pair<gemmi::Atom*, gemmi::Residue*>> sel;
    for (const gemmi::Chain& chain : model.chains) {
        if (chain.name == chain_name) {
            for (const gemmi::Residue& res : chain.residues) {
                for (const gemmi::Atom& atom : res.atoms) {
                    sel.push_back({const_cast<gemmi::Atom*>(&atom),
                                   const_cast<gemmi::Residue*>(&res)});
                }
            }
        }
    }
    return sel;
}

TEST(OracleTest, merge_atom_selections) {
    // Case 1: Downstream short fragment merge — Chain A has res 1-5, Chain B has res 5-7
    // B is shorter (3 vs 5 residues), overlap at residue 5 (atoms within 0.02Å)
    {
        std::ofstream pdb1("/tmp/fragment_test.pdb");
        char buf[80];
        int serial = 1;
        for (int i = 1; i <= 5; i++) {
            float x = 1.0f + i;
            sprintf(buf, "ATOM  %5d  CA  ALA A%4d    %8.3f%8.3f%8.3f  1.00 20.00\n", serial++, i, x, 2.0f, 3.0f); pdb1 << buf;
            sprintf(buf, "ATOM  %5d  N   ALA A%4d    %8.3f%8.3f%8.3f  1.00 20.00\n", serial++, i, x, 2.5f, 3.5f); pdb1 << buf;
            sprintf(buf, "ATOM  %5d  C   ALA A%4d    %8.3f%8.3f%8.3f  1.00 20.00\n", serial++, i, x+1.0f, 2.0f, 3.0f); pdb1 << buf;
            sprintf(buf, "ATOM  %5d  O   ALA A%4d    %8.3f%8.3f%8.3f  1.00 20.00\n", serial++, i, x+1.5f, 2.5f, 3.5f); pdb1 << buf;
        }
        for (int i = 5; i <= 7; i++) {
            float x = 1.0f + i;
            sprintf(buf, "ATOM  %5d  CA  ALA B%4d    %8.3f%8.3f%8.3f  1.00 20.00\n", serial++, i, x+0.01f, 2.01f, 3.01f); pdb1 << buf;
            sprintf(buf, "ATOM  %5d  N   ALA B%4d    %8.3f%8.3f%8.3f  1.00 20.00\n", serial++, i, x+0.01f, 2.51f, 3.51f); pdb1 << buf;
            sprintf(buf, "ATOM  %5d  C   ALA B%4d    %8.3f%8.3f%8.3f  1.00 20.00\n", serial++, i, x+1.01f, 2.01f, 3.01f); pdb1 << buf;
            sprintf(buf, "ATOM  %5d  O   ALA B%4d    %8.3f%8.3f%8.3f  1.00 20.00\n", serial++, i, x+1.51f, 2.51f, 3.51f); pdb1 << buf;
        }
        pdb1.close();

        gemmi::Structure st = gemmi::read_pdb_file("/tmp/fragment_test.pdb");

        EXPECT_FALSE(st.models.empty());
        gemmi::Model& model = st.models[0];

        int n_chains_before = static_cast<int>(model.chains.size());
        EXPECT_EQ(n_chains_before, 2);

        auto sel1 = build_selection(model, "A");
        auto sel2 = build_selection(model, "B");

        EXPECT_EQ(static_cast<int>(sel1.size()), 20);
        EXPECT_EQ(static_cast<int>(sel2.size()), 12);

        coot::delete_a_chain_t result = coot::merge_atom_selections_gemmi(st, sel1, sel2);

        EXPECT_TRUE(result.chains_were_mergeable);
        EXPECT_FALSE(result.short_fragment_is_in_first_selection);
        EXPECT_TRUE(result.short_fragment_is_upstream_fragment);
        EXPECT_EQ(static_cast<int>(result.delete_type), 0);
        EXPECT_EQ(static_cast<int>(model.chains.size()), 2);
    }

    // Case 2: Non-overlapping — protein chain A vs HOH chain B from standard PDB
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        EXPECT_FALSE(st.models.empty());
        gemmi::Model& model = st.models[0];

        auto sel1 = build_selection(model, "A");
        auto sel2 = build_selection(model, "B");

        coot::delete_a_chain_t result = coot::merge_atom_selections_gemmi(st, sel1, sel2);

        EXPECT_FALSE(result.chains_were_mergeable);
        EXPECT_FALSE(result.short_fragment_is_in_first_selection);
        EXPECT_FALSE(result.short_fragment_is_upstream_fragment);
    }

    // Case 3: Same atoms but too far apart (>0.8 Å) — should not merge
    {
        std::ofstream pdb3("/tmp/no_overlap_test.pdb");
        char buf[80];
        int serial = 1;
        for (int i = 1; i <= 5; i++) {
            float x = 1.0f + i;
            sprintf(buf, "ATOM  %5d  CA  ALA A%4d    %8.3f%8.3f%8.3f  1.00 20.00\n", serial++, i, x, 2.0f, 3.0f); pdb3 << buf;
            sprintf(buf, "ATOM  %5d  N   ALA A%4d    %8.3f%8.3f%8.3f  1.00 20.00\n", serial++, i, x, 2.5f, 3.5f); pdb3 << buf;
            sprintf(buf, "ATOM  %5d  C   ALA A%4d    %8.3f%8.3f%8.3f  1.00 20.00\n", serial++, i, x+1.0f, 2.0f, 3.0f); pdb3 << buf;
            sprintf(buf, "ATOM  %5d  O   ALA A%4d    %8.3f%8.3f%8.3f  1.00 20.00\n", serial++, i, x+1.5f, 2.5f, 3.5f); pdb3 << buf;
        }
        for (int i = 5; i <= 7; i++) {
            float x = 21.0f + i;
            sprintf(buf, "ATOM  %5d  CA  ALA B%4d    %8.3f%8.3f%8.3f  1.00 20.00\n", serial++, i, x, 2.0f, 3.0f); pdb3 << buf;
            sprintf(buf, "ATOM  %5d  N   ALA B%4d    %8.3f%8.3f%8.3f  1.00 20.00\n", serial++, i, x, 2.5f, 3.5f); pdb3 << buf;
            sprintf(buf, "ATOM  %5d  C   ALA B%4d    %8.3f%8.3f%8.3f  1.00 20.00\n", serial++, i, x+1.0f, 2.0f, 3.0f); pdb3 << buf;
            sprintf(buf, "ATOM  %5d  O   ALA B%4d    %8.3f%8.3f%8.3f  1.00 20.00\n", serial++, i, x+1.5f, 2.5f, 3.5f); pdb3 << buf;
        }
        pdb3.close();

        gemmi::Structure st = gemmi::read_pdb_file("/tmp/no_overlap_test.pdb");

        EXPECT_FALSE(st.models.empty());
        gemmi::Model& model = st.models[0];

        auto sel1 = build_selection(model, "A");
        auto sel2 = build_selection(model, "B");

        coot::delete_a_chain_t result = coot::merge_atom_selections_gemmi(st, sel1, sel2);

        EXPECT_FALSE(result.chains_were_mergeable);
        EXPECT_FALSE(result.short_fragment_is_in_first_selection);
        EXPECT_FALSE(result.short_fragment_is_upstream_fragment);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}