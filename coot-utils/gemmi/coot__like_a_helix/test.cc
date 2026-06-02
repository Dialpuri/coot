#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <map>
#include "function.hh"

TEST(OracleTest, like_a_helix) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    ASSERT_FALSE(st.models.empty());

    // Find chain A
    gemmi::Model* model = &st.models[0];
    gemmi::Chain* chain_a = nullptr;
    for (gemmi::Chain& ch : model->chains) {
        if (ch.name == "A") {
            chain_a = &ch;
            break;
        }
    }
    ASSERT_NE(chain_a, nullptr);

    // Build a map from residue pointer to chain name for assertions
    std::map<const gemmi::Residue*, std::string> residue_to_chain;
    for (gemmi::Chain& ch : model->chains) {
        for (gemmi::Residue& res : ch.residues) {
            residue_to_chain[&res] = ch.name;
        }
    }

    // Case 1: select all residues in chain A (full chain)
    {
        std::vector<gemmi::Residue*> selected_residues;
        for (gemmi::Residue& res : chain_a->residues) {
            selected_residues.push_back(&res);
        }

        std::vector<gemmi::Residue*> helical = coot::like_a_helix_gemmi(selected_residues);

        EXPECT_EQ(helical.size(), 54u);

        // First 10 residues
        EXPECT_NE(helical[0], nullptr);
        EXPECT_EQ(residue_to_chain[helical[0]], std::string("A"));
        EXPECT_EQ(helical[0]->seqid.num.value, 47);
        EXPECT_EQ(helical[0]->name, std::string("THR"));

        EXPECT_NE(helical[1], nullptr);
        EXPECT_EQ(residue_to_chain[helical[1]], std::string("A"));
        EXPECT_EQ(helical[1]->seqid.num.value, 48);
        EXPECT_EQ(helical[1]->name, std::string("ALA"));

        EXPECT_NE(helical[2], nullptr);
        EXPECT_EQ(residue_to_chain[helical[2]], std::string("A"));
        EXPECT_EQ(helical[2]->seqid.num.value, 49);
        EXPECT_EQ(helical[2]->name, std::string("ILE"));

        EXPECT_NE(helical[3], nullptr);
        EXPECT_EQ(residue_to_chain[helical[3]], std::string("A"));
        EXPECT_EQ(helical[3]->seqid.num.value, 50);
        EXPECT_EQ(helical[3]->name, std::string("ARG"));

        EXPECT_NE(helical[4], nullptr);
        EXPECT_EQ(residue_to_chain[helical[4]], std::string("A"));
        EXPECT_EQ(helical[4]->seqid.num.value, 51);
        EXPECT_EQ(helical[4]->name, std::string("GLU"));

        EXPECT_NE(helical[5], nullptr);
        EXPECT_EQ(residue_to_chain[helical[5]], std::string("A"));
        EXPECT_EQ(helical[5]->seqid.num.value, 52);
        EXPECT_EQ(helical[5]->name, std::string("ILE"));

        EXPECT_NE(helical[6], nullptr);
        EXPECT_EQ(residue_to_chain[helical[6]], std::string("A"));
        EXPECT_EQ(helical[6]->seqid.num.value, 53);
        EXPECT_EQ(helical[6]->name, std::string("SER"));

        EXPECT_NE(helical[7], nullptr);
        EXPECT_EQ(residue_to_chain[helical[7]], std::string("A"));
        EXPECT_EQ(helical[7]->seqid.num.value, 54);
        EXPECT_EQ(helical[7]->name, std::string("LEU"));

        EXPECT_NE(helical[8], nullptr);
        EXPECT_EQ(residue_to_chain[helical[8]], std::string("A"));
        EXPECT_EQ(helical[8]->seqid.num.value, 88);
        EXPECT_EQ(helical[8]->name, std::string("LYS"));

        EXPECT_NE(helical[9], nullptr);
        EXPECT_EQ(residue_to_chain[helical[9]], std::string("A"));
        EXPECT_EQ(helical[9]->seqid.num.value, 89);
        EXPECT_EQ(helical[9]->name, std::string("LYS"));
    }

    // Case 2: select a very short region (too few residues for helix detection)
    // Residues 2-3 in chain A (only 2 residues)
    {
        std::vector<gemmi::Residue*> selected_residues;
        for (gemmi::Residue& res : chain_a->residues) {
            if (res.seqid.num.value >= 2 && res.seqid.num.value <= 3) {
                selected_residues.push_back(&res);
            }
        }

        std::vector<gemmi::Residue*> helical = coot::like_a_helix_gemmi(selected_residues);

        EXPECT_EQ(helical.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}