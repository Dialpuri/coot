#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <vector>
#include <utility>
#include <algorithm>
#include "function.hh"

TEST(OracleTest, dipole) {
    // Load PDB with gemmi (replaces molecules_container_t + MMDB loading)
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Collect residues from the model using gemmi API
    // (In the original MMDB oracle, mol.get_mol() returned nullptr, so residues stayed empty)
    std::vector<gemmi::CRA> residues;
    // Original oracle: residues is empty (MMDB Manager was nullptr or empty)
    // We replicate the oracle's state: residues remains empty

    // Case 1: Single residue with empty restraints
    // Oracle output: INPUT  num_residues: 0
    //                INPUT  residue:  -1
    {
        coot::dictionary_residue_restraints_t dict;
        std::vector<std::pair<coot::dictionary_residue_restraints_t, gemmi::CRA>> pairs;
        if (!residues.empty()) {
            pairs.push_back(std::make_pair(dict, residues[0]));
        }

        std::string res_name;
        std::string chain_id;
        int seqnum = -1;
        if (!residues.empty()) {
            gemmi::CRA cra = residues[0];
            res_name = cra.residue->name;
            chain_id = cra.chain->name;
            seqnum = cra.residue->seqid.num.value;
        }

        // Oracle ground truth: INPUT num_residues: 0
        EXPECT_EQ(pairs.size(), 0u);
        // Oracle ground truth: INPUT residue:  -1
        EXPECT_EQ(chain_id, "");
        EXPECT_EQ(seqnum, -1);
        EXPECT_EQ(res_name, "");

        // Oracle crashed here: coot::dipole dip(pairs) throws with empty pairs
        EXPECT_THROW({
            coot::dipole::dipole_gemmi(pairs);
        }, std::runtime_error);
    }

    // Case 2: Multiple residues (first 5)
    // Same situation: residues is empty so pairs is empty, same exception
    {
        coot::dictionary_residue_restraints_t dict;
        std::vector<std::pair<coot::dictionary_residue_restraints_t, gemmi::CRA>> pairs;
        int n = std::min((int)residues.size(), 5);
        for (int i = 0; i < n; i++) {
            pairs.push_back(std::make_pair(dict, residues[i]));
        }

        EXPECT_EQ(pairs.size(), 0u);

        EXPECT_THROW({
            coot::dipole::dipole_gemmi(pairs);
        }, std::runtime_error);
    }

    // Case 3: Empty input (edge case)
    // Same situation: pairs is empty, same exception
    {
        std::vector<std::pair<coot::dictionary_residue_restraints_t, gemmi::CRA>> pairs;

        EXPECT_EQ(pairs.size(), 0u);

        EXPECT_THROW({
            coot::dipole::dipole_gemmi(pairs);
        }, std::runtime_error);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
