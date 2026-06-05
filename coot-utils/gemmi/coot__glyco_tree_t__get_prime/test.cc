#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "coot-utils/tree.hh"
#include <algorithm>
#include "function.hh"

static gemmi::Structure* load_structure() {
    gemmi::Structure* st = nullptr;
    std::string path = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";
    try { st = new gemmi::Structure(gemmi::read_pdb_file(path)); }
    catch (...) { st = nullptr; }
    if (!st) {
        path = "/lmb/home/ds5/test-data/example.pdb";
        try { st = new gemmi::Structure(gemmi::read_pdb_file(path)); }
        catch (...) { st = nullptr; }
    }
    if (!st) {
        path = "/tmp/example.pdb";
        try { st = new gemmi::Structure(gemmi::read_pdb_file(path)); }
        catch (...) { st = nullptr; }
    }
    return st;
}

TEST(OracleTest, glyco_tree_t_get_prime) {
    // case: ALPHA1-3 child -> NON_PRIME
    {
        gemmi::Structure* st = load_structure();
        ASSERT_NE(st, nullptr);

        EXPECT_FALSE(st->models.empty());
        gemmi::Model& model = st->models[0];

        auto chain_it = std::find_if(model.chains.begin(), model.chains.end(),
            [](const gemmi::Chain& c) { return c.name == "A"; });
        ASSERT_NE(chain_it, model.chains.end());
        gemmi::Chain& chain = *chain_it;

        ASSERT_GE(chain.residues.size(), 6u);

        gemmi::Residue* res_bma = &chain.residues[0];
        gemmi::Residue* res_alpha13 = &chain.residues[1];
        gemmi::Residue* res_alpha16 = &chain.residues[2];

        // Build tree manually - set root using tree constructor
        coot::linked_residue_t_gemmi bma_data(res_bma, "BMA", "");
        tree<coot::linked_residue_t_gemmi> glyco_tree(bma_data);
        tree<coot::linked_residue_t_gemmi>::iterator it_root = glyco_tree.begin();

        coot::linked_residue_t_gemmi a13_data(res_alpha13, "MAN", "ALPHA1-3");
        glyco_tree.append_child(it_root, a13_data);

        coot::linked_residue_t_gemmi a16_data(res_alpha16, "MAN", "ALPHA1-6");
        glyco_tree.append_child(it_root, a16_data);

        int result = coot::get_prime_gemmi(glyco_tree, res_alpha13);
        EXPECT_EQ(result, 2); // NON_PRIME

        delete st;
    }

    // case: ALPHA1-6 child -> PRIME
    {
        gemmi::Structure* st = load_structure();
        ASSERT_NE(st, nullptr);

        EXPECT_FALSE(st->models.empty());
        gemmi::Model& model = st->models[0];

        auto chain_it = std::find_if(model.chains.begin(), model.chains.end(),
            [](const gemmi::Chain& c) { return c.name == "A"; });
        ASSERT_NE(chain_it, model.chains.end());
        gemmi::Chain& chain = *chain_it;

        ASSERT_GE(chain.residues.size(), 6u);

        gemmi::Residue* res_bma = &chain.residues[0];
        gemmi::Residue* res_alpha13 = &chain.residues[1];
        gemmi::Residue* res_alpha16 = &chain.residues[2];

        // Build tree manually - set root using tree constructor
        coot::linked_residue_t_gemmi bma_data(res_bma, "BMA", "");
        tree<coot::linked_residue_t_gemmi> glyco_tree(bma_data);
        tree<coot::linked_residue_t_gemmi>::iterator it_root = glyco_tree.begin();

        coot::linked_residue_t_gemmi a13_data(res_alpha13, "MAN", "ALPHA1-3");
        glyco_tree.append_child(it_root, a13_data);

        coot::linked_residue_t_gemmi a16_data(res_alpha16, "MAN", "ALPHA1-6");
        glyco_tree.append_child(it_root, a16_data);

        int result = coot::get_prime_gemmi(glyco_tree, res_alpha16);
        EXPECT_EQ(result, 1); // PRIME

        delete st;
    }

    // case: residue not in tree -> UNSET
    {
        gemmi::Structure* st = load_structure();
        ASSERT_NE(st, nullptr);

        EXPECT_FALSE(st->models.empty());
        gemmi::Model& model = st->models[0];

        auto chain_it = std::find_if(model.chains.begin(), model.chains.end(),
            [](const gemmi::Chain& c) { return c.name == "A"; });
        ASSERT_NE(chain_it, model.chains.end());
        gemmi::Chain& chain = *chain_it;

        ASSERT_GE(chain.residues.size(), 6u);

        gemmi::Residue* res_bma = &chain.residues[0];
        gemmi::Residue* res_alpha13 = &chain.residues[1];
        gemmi::Residue* res_alpha16 = &chain.residues[2];

        // Build tree manually - set root using tree constructor
        coot::linked_residue_t_gemmi bma_data(res_bma, "BMA", "");
        tree<coot::linked_residue_t_gemmi> glyco_tree(bma_data);
        tree<coot::linked_residue_t_gemmi>::iterator it_root = glyco_tree.begin();

        coot::linked_residue_t_gemmi a13_data(res_alpha13, "MAN", "ALPHA1-3");
        glyco_tree.append_child(it_root, a13_data);

        coot::linked_residue_t_gemmi a16_data(res_alpha16, "MAN", "ALPHA1-6");
        glyco_tree.append_child(it_root, a16_data);

        gemmi::Residue* res_not_in_tree = &chain.residues[5];
        ASSERT_NE(res_not_in_tree, nullptr);

        int result = coot::get_prime_gemmi(glyco_tree, res_not_in_tree);
        EXPECT_EQ(result, 0); // UNSET

        delete st;
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}