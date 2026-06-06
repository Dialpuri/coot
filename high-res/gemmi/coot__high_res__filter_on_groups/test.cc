#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: collect all atom CRAs from a model
static std::vector<gemmi::CRA> collect_all_atom_cras(gemmi::Model& model) {
    std::vector<gemmi::CRA> cras;
    for (auto& chain : model.chains) {
        for (auto& res : chain.residues) {
            for (auto& atom : res.atoms) {
                cras.push_back(gemmi::CRA{&chain, &res, &atom});
            }
        }
    }
    return cras;
}

TEST(OracleTest, filter_on_groups) {
    // Load PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    ASSERT_GT(st.models.size(), 0u) << "No models loaded";

    gemmi::Model& model = st.models[0];

    // Create atom selection for all atoms as CRAs
    std::vector<gemmi::CRA> atom_selection = collect_all_atom_cras(model);

    ASSERT_GT(atom_selection.size(), 0u) << "No atoms selected";

    // Case 1: Simple group with 3 atoms (first 3 atoms as one group)
    {
        std::vector<std::vector<int>> groups;
        groups.push_back({0, 1, 2});

        coot::minimol::molecule result = coot::high_res::filter_on_groups_gemmi(
            groups, atom_selection);

        EXPECT_EQ(result.fragments.size(), 1u);
        if (!result.fragments.empty()) {
            EXPECT_EQ(result.fragments[0].residues.size(), 2u);
            if (!result.fragments[0].residues.empty()) {
                auto &res = result.fragments[0].residues[0];
                EXPECT_EQ(res.name, "");
                EXPECT_EQ(res.seqnum, 1);
                EXPECT_EQ(res.atoms.size(), 0u);
            }
        }
    }

    // Case 2: Two groups
    {
        std::vector<std::vector<int>> groups;
        groups.push_back({0, 1});
        groups.push_back({3, 4, 5});

        coot::minimol::molecule result = coot::high_res::filter_on_groups_gemmi(
            groups, atom_selection);

        EXPECT_EQ(result.fragments.size(), 1u);
        if (!result.fragments.empty()) {
            EXPECT_EQ(result.fragments[0].residues.size(), 2u);
            if (!result.fragments[0].residues.empty()) {
                auto &res = result.fragments[0].residues[0];
                EXPECT_EQ(res.name, "");
                EXPECT_EQ(res.seqnum, 1);
                EXPECT_EQ(res.atoms.size(), 0u);
            }
        }
    }

    // Case 3: Empty groups vector (edge case)
    {
        std::vector<std::vector<int>> groups;

        coot::minimol::molecule result = coot::high_res::filter_on_groups_gemmi(
            groups, atom_selection);

        EXPECT_EQ(result.fragments.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}