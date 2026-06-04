#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, atoms_match_dictionary) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    ASSERT_FALSE(st.models.empty());

    gemmi::Model& model = st.models[0];
    ASSERT_FALSE(model.chains.empty());

    // Get chain A
    gemmi::Chain* chain = nullptr;
    for (gemmi::Chain& ch : model.chains) {
        if (ch.name == "A") {
            chain = &ch;
            break;
        }
    }
    ASSERT_NE(chain, nullptr);

    int n_residues = static_cast<int>(chain->residues.size());

    // Case 1: valid protein residues (first 3 residues of chain A)
    {
        std::vector<gemmi::Residue*> residues;
        int num_test = std::min(n_residues, 3);
        for (int i = 0; i < num_test; i++) {
            residues.push_back(&chain->residues[i]);
        }

        bool check_hydrogens_too_flag = false;
        bool apply_bond_distance_check = false;

        auto result = coot::atoms_match_dictionary_gemmi(residues, check_hydrogens_too_flag, apply_bond_distance_check);

        EXPECT_TRUE(result.first);
        EXPECT_EQ(result.second.size(), 0u);
    }

    // Case 2: empty residue vector
    {
        std::vector<gemmi::Residue*> residues;

        auto result = coot::atoms_match_dictionary_gemmi(residues, false, false);

        EXPECT_TRUE(result.first);
        EXPECT_EQ(result.second.size(), 0u);
    }

    // Case 3: single standard residue with hydrogen check enabled
    {
        std::vector<gemmi::Residue*> residues;
        residues.push_back(&chain->residues[0]);

        bool check_hydrogens_too_flag = true;
        bool apply_bond_distance_check = true;

        auto result = coot::atoms_match_dictionary_gemmi(residues, check_hydrogens_too_flag, apply_bond_distance_check);

        EXPECT_TRUE(result.first);
        EXPECT_EQ(result.second.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}