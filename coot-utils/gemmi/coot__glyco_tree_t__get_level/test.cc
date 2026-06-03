#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

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

TEST(OracleTest, glyco_tree_t_get_level_gemmi) {
    {
        // Case: root, child, grandchild, and not_in_tree levels
        gemmi::Structure st =
            gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        ASSERT_FALSE(st.models.empty());
        gemmi::Model& model = st.models[0];

        // Get residues from chain A: res50, res51, res52
        gemmi::Residue* res1 = find_residue_in_model(model, "A", 50);
        gemmi::Residue* res2 = find_residue_in_model(model, "A", 51);
        gemmi::Residue* res3 = find_residue_in_model(model, "A", 52);

        ASSERT_NE(res1, nullptr);
        ASSERT_NE(res2, nullptr);
        ASSERT_NE(res3, nullptr);

        // Build tree manually using gemmi-linked_residue_t: res1 as root, res2 as child
        // of res1, res3 as child of res2
        coot::linked_residue_t_gemmi lr1;
        lr1.residue = res1;
        lr1.link_type = "";

        coot::linked_residue_t_gemmi lr2;
        lr2.residue = res2;
        lr2.link_type = "beta1-4";

        coot::linked_residue_t_gemmi lr3;
        lr3.residue = res3;
        lr3.link_type = "beta1-4";

        // Build the tree: root with res1, child with res2, grandchild with res3
        tree<coot::linked_residue_t_gemmi> my_tree;
        my_tree.insert(my_tree.end(), lr1);
        tree<coot::linked_residue_t_gemmi>::iterator it1 = my_tree.begin();
        my_tree.append_child(it1, lr2);
        tree<coot::linked_residue_t_gemmi>::iterator it2 = it1;
        ++it2;
        my_tree.append_child(it2, lr3);

        // Case 1: root node (no parent) -> level should be 0
        int level1 = coot::get_level_gemmi(my_tree, res1);
        EXPECT_EQ(level1, 0);

        // Case 2: child of root -> level should be 1
        int level2 = coot::get_level_gemmi(my_tree, res2);
        EXPECT_EQ(level2, 1);

        // Case 3: grandchild -> level should be 2
        int level3 = coot::get_level_gemmi(my_tree, res3);
        EXPECT_EQ(level3, 2);

        // Edge case: residue not in tree -> level should be -1
        gemmi::Residue* res_not = find_residue_in_model(model, "A", 53);
        if (res_not) {
            int level_not = coot::get_level_gemmi(my_tree, res_not);
            EXPECT_EQ(level_not, -1);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}