#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include "coot-utils/tree.hh"

static gemmi::Residue* find_residue_in_model(gemmi::Model& model,
                                             const std::string& chain_id,
                                             int seq_num) {
    for (auto& chain : model.chains) {
        if (chain.name == chain_id) {
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == seq_num) {
                    return &res;
                }
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, glyco_tree_t_get_id) {
    gemmi::Structure st =
        gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    ASSERT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];

    // --- Case 1: a valid residue (protein residue, no glycan linkage) ---
    // Use a residue that exists in the tree but has no parent (root only)
    // so get_id returns the default-empty result matching the oracle
    {
        gemmi::Residue* first_res = find_residue_in_model(model, "A", 50);
        ASSERT_NE(first_res, nullptr);

        EXPECT_EQ("A", "A");
        EXPECT_EQ(first_res->seqid.num.value, 50);

        // Build tree with single root node (no parent → empty id fields)
        coot::linked_residue_t_gemmi lr1;
        lr1.residue = first_res;
        lr1.link_type = "";

        tree<coot::linked_residue_t_gemmi> glyco_tree;
        glyco_tree.insert(glyco_tree.end(), lr1);

        coot::residue_id_t_gemmi id = coot::get_id_gemmi(glyco_tree, first_res);

        EXPECT_EQ(id.level, 0);
        EXPECT_EQ(id.res_type, "");
        EXPECT_EQ(id.link_type, "");
        EXPECT_EQ(id.parent_res_type, "");
    }

    // --- Case 2: null residue pointer ---
    {
        tree<coot::linked_residue_t_gemmi> empty_tree;

        coot::residue_id_t_gemmi id = coot::get_id_gemmi(empty_tree, nullptr);

        EXPECT_EQ(id.level, 0);
        EXPECT_EQ(id.res_type, "");
    }

    // --- Case 3: residue WITH a parent in the tree (non-trivial path) ---
    {
        gemmi::Residue* res1 = find_residue_in_model(model, "A", 50);
        gemmi::Residue* res2 = find_residue_in_model(model, "A", 51);

        ASSERT_NE(res1, nullptr);
        ASSERT_NE(res2, nullptr);

        // Build tree: res1 as root, res2 as child with linkage
        coot::linked_residue_t_gemmi lr1;
        lr1.residue = res1;
        lr1.link_type = "";

        coot::linked_residue_t_gemmi lr2;
        lr2.residue = res2;
        lr2.link_type = "beta1-4";

        tree<coot::linked_residue_t_gemmi> glyco_tree;
        glyco_tree.insert(glyco_tree.end(), lr1);
        tree<coot::linked_residue_t_gemmi>::iterator it1 = glyco_tree.begin();
        glyco_tree.append_child(it1, lr2);

        // Get id for the child (res2) — should have non-empty parent info
        coot::residue_id_t_gemmi id = coot::get_id_gemmi(glyco_tree, res2);

        EXPECT_EQ(id.level, 1);
        EXPECT_EQ(id.res_type, res2->name);
        EXPECT_EQ(id.link_type, "beta1-4");
        EXPECT_EQ(id.parent_res_type, res1->name);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}