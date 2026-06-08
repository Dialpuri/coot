#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_asc_index_new) {
    // Read the test PDB with gemmi
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Find chain A, residue 10 — this is what the original test uses
    gemmi::Residue* target_res = nullptr;
    for (gemmi::Model& model : st.models) {
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name != "A") continue;
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == 10) {
                    target_res = &res;
                    break;
                }
            }
        }
    }
    ASSERT_NE(target_res, nullptr);

    std::vector<std::pair<bool, gemmi::Residue*>> residues;
    residues.push_back(std::make_pair(false, target_res));

    coot::restraints_container_t_gemmi restraints(residues, &st);

    // Case 1: valid atom — CA in residue 10, chain A
    {
        int index = coot::get_asc_index_new_gemmi(
            restraints, &st, "CA", "", 10, "", "A");
        EXPECT_EQ(index, 1);
    }

    // Case 2: invalid atom name — non-existent atom
    {
        int index = coot::get_asc_index_new_gemmi(
            restraints, &st, "ZZZ", "", 10, "", "A");
        EXPECT_EQ(index, -1);
    }

    // Case 3: non-existent residue
    {
        int index = coot::get_asc_index_new_gemmi(
            restraints, &st, "CA", "", 9999, "", "A");
        EXPECT_EQ(index, -1);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
