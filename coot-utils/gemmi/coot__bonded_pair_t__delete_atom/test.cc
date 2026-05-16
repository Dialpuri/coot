#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper to count atoms on a gemmi residue.
static int count_atoms(gemmi::Residue *res) {
    return static_cast<int>(res->atoms.size());
}

TEST(OracleTest, bonded_pair_delete_atom) {
    // Load the PDB structure.
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Find chain A, residue 10.
    gemmi::Residue *res = nullptr;
    for (auto &model : st.models) {
        for (auto &chain : model.chains) {
            if (chain.name == "A") {
                for (auto &r : chain.residues) {
                    if (r.seqid.num.value == 10 && r.seqid.icode == ' ') {
                        res = &r;
                        break;
                    }
                }
            }
            if (res) break;
        }
        if (res) break;
    }
    ASSERT_NE(res, nullptr);

    std::string chain_id = "A";
    std::string res_name = res->name;
    int seq_num = res->seqid.num.value;
    std::string residue_info = chain_id + " " + std::to_string(seq_num) + " " + res_name;

    // --- Case 1: delete an existing atom ("CB" from ILE 10) ---
    {
        int n_before = count_atoms(res);
        std::string atom_name = "CB";

        EXPECT_EQ(residue_info, "A 10 ILE");
        EXPECT_EQ(atom_name, "CB");
        EXPECT_EQ(n_before, 8);

        coot::bonded_pair_t bp;
        bp.delete_atom_gemmi(res, atom_name);

        int n_after = count_atoms(res);
        EXPECT_EQ(n_after, 7);
    }

    // --- Case 2: delete a non-existent atom (no-op) ---
    {
        int n_before = count_atoms(res);
        std::string atom_name = "ZZZ";

        EXPECT_EQ(n_before, 7);

        coot::bonded_pair_t bp;
        bp.delete_atom_gemmi(res, atom_name);

        int n_after = count_atoms(res);
        EXPECT_EQ(n_after, 7);
    }

    // --- Case 3: delete "CG1" (another existing atom) ---
    {
        int n_before = count_atoms(res);
        std::string atom_name = "CG1";

        EXPECT_EQ(n_before, 7);

        coot::bonded_pair_t bp;
        bp.delete_atom_gemmi(res, atom_name);

        int n_after = count_atoms(res);
        EXPECT_EQ(n_after, 6);
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}