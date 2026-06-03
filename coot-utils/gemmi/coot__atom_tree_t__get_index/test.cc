#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, atom_tree_t_get_index) {
    // Load the PDB file with gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Get model 0 (gemmi is 0-indexed)
    gemmi::Model& model = st.models[0];

    // Find chain A
    gemmi::Chain* chain_a = nullptr;
    for (gemmi::Chain& ch : model.chains) {
        if (ch.name == "A") {
            chain_a = &ch;
            break;
        }
    }
    ASSERT_NE(chain_a, nullptr);

    // Find residue //A/1 (MET, seqid.num.value == 1)
    gemmi::Residue* res = nullptr;
    for (gemmi::Residue& r : chain_a->residues) {
        if (r.seqid.num.value == 1) {
            res = &r;
            break;
        }
    }
    ASSERT_NE(res, nullptr);

    // Get atoms from residue
    ASSERT_GE(res->atoms.size(), 2u);

    // Build atom_tree_t using contact_indices constructor (gemmi variant, sets residue_g internally)
    std::vector<std::vector<int>> empty_contacts;
    int base_atom_index = 0;
    std::string altconf;

    coot::atom_tree_t tree(empty_contacts, base_atom_index, res, altconf);

    // Case 1: valid atom in the residue (first atom, N)
    {
        gemmi::Atom* first_atom = &(res->atoms[0]);
        coot::map_index_t idx = tree.get_index_gemmi(first_atom);
        EXPECT_EQ(idx.index(), 0);
    }

    // Case 2: second atom in the same residue (CA)
    {
        gemmi::Atom* second_atom = &(res->atoms[1]);
        coot::map_index_t idx2 = tree.get_index_gemmi(second_atom);
        EXPECT_EQ(idx2.index(), 1);
    }

    // Case 3: atom NOT in the residue (from a different residue //A/2)
    {
        // Find residue //A/2
        gemmi::Residue* res2 = nullptr;
        for (gemmi::Residue& r : chain_a->residues) {
            if (r.seqid.num.value == 2) {
                res2 = &r;
                break;
            }
        }
        if (res2 && !res2->atoms.empty()) {
            gemmi::Atom* other_atom = &(res2->atoms[0]);
            coot::map_index_t idx3 = tree.get_index_gemmi(other_atom);
            EXPECT_EQ(idx3.index(), -1);
        }
    }

    // Case 4: null atom
    {
        coot::map_index_t idx4 = tree.get_index_gemmi(nullptr);
        EXPECT_EQ(idx4.index(), -1);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}