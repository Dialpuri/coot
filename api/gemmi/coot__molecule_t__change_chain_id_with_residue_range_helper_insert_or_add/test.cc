#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, change_chain_id_with_residue_range_helper_insert_or_add) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model& model = st.models[0];

    // Get chain A
    gemmi::Chain* chain_a = nullptr;
    for (gemmi::Chain& c : model.chains) {
        if (c.name == "A") {
            chain_a = &c;
            break;
        }
    }
    ASSERT_NE(chain_a, nullptr);

    // Get residue from chain A to copy (used across all cases)
    ASSERT_FALSE(chain_a->residues.empty());
    gemmi::Residue orig_res = chain_a->residues[0];

    // --- Case 1: insert_middle ---
    {
        int copy_seqnum = 5;
        gemmi::Residue res_copy = orig_res;
        res_copy.seqid.num.value = copy_seqnum;

        int residues_before = static_cast<int>(chain_a->residues.size());
        EXPECT_EQ(residues_before, 267);

        coot::molecule_t mol;
        mol.change_chain_id_with_residue_range_helper_insert_or_add_gemmi(*chain_a, res_copy);

        int residues_after = static_cast<int>(chain_a->residues.size());
        EXPECT_EQ(residues_after, 268);
    }

    // --- Case 2: add_end ---
    {
        int copy_seqnum = 9999;
        gemmi::Residue res_copy2 = orig_res;
        res_copy2.seqid.num.value = copy_seqnum;

        int residues_before = static_cast<int>(chain_a->residues.size());
        EXPECT_EQ(residues_before, 268);

        coot::molecule_t mol;
        mol.change_chain_id_with_residue_range_helper_insert_or_add_gemmi(*chain_a, res_copy2);

        int residues_after = static_cast<int>(chain_a->residues.size());
        EXPECT_EQ(residues_after, 269);
    }

    // --- Case 3: empty_chain ---
    {
        gemmi::Chain new_chain("Z");
        model.chains.push_back(new_chain);

        gemmi::Residue res_copy3 = orig_res;

        gemmi::Chain* new_chain_ptr = &model.chains.back();
        int residues_before = static_cast<int>(new_chain_ptr->residues.size());
        EXPECT_EQ(residues_before, 0);

        coot::molecule_t mol;
        mol.change_chain_id_with_residue_range_helper_insert_or_add_gemmi(*new_chain_ptr, res_copy3);

        int residues_after = static_cast<int>(new_chain_ptr->residues.size());
        EXPECT_EQ(residues_after, 1);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}