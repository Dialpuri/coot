#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_mainchain_atoms) {
    // Setup
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model& model = st.models[0];

    // Get chain A
    const gemmi::Chain* chain_A = nullptr;
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == "A") {
            chain_A = &chain;
            break;
        }
    }
    ASSERT_NE(chain_A, nullptr);

    // Get residues from chain A
    ASSERT_GE(chain_A->residues.size(), 2u);

    const gemmi::Residue& res_1 = chain_A->residues[0];
    const gemmi::Residue& res_2 = chain_A->residues[1];

    // Verify inputs match oracle
    EXPECT_EQ(std::string(chain_A->name), "A");
    EXPECT_EQ(res_1.seqid.num.value, 1);
    EXPECT_EQ(std::string(res_1.name), "MET");
    EXPECT_EQ(std::string(chain_A->name), "A");
    EXPECT_EQ(res_2.seqid.num.value, 2);
    EXPECT_EQ(std::string(res_2.name), "GLU");

    // Case 1: two consecutive residues from same chain
    {
        std::vector<const gemmi::Atom*> mainchain_atoms = coot::crankshaft::get_mainchain_atoms_gemmi(res_1, res_2);
        EXPECT_EQ(mainchain_atoms.size(), 4u);

        EXPECT_EQ(std::string(mainchain_atoms[0]->name), "N");
        EXPECT_EQ(res_1.seqid.num.value, 1);
        EXPECT_EQ(std::string(res_1.name), "MET");

        EXPECT_EQ(std::string(mainchain_atoms[1]->name), "CA");
        EXPECT_EQ(res_1.seqid.num.value, 1);
        EXPECT_EQ(std::string(res_1.name), "MET");

        EXPECT_EQ(std::string(mainchain_atoms[2]->name), "C");
        EXPECT_EQ(res_1.seqid.num.value, 1);
        EXPECT_EQ(std::string(res_1.name), "MET");

        EXPECT_EQ(std::string(mainchain_atoms[3]->name), "O");
        EXPECT_EQ(res_1.seqid.num.value, 1);
        EXPECT_EQ(std::string(res_1.name), "MET");
    }

    // Case 2: same residue twice
    {
        std::vector<const gemmi::Atom*> mc_same = coot::crankshaft::get_mainchain_atoms_gemmi(res_1, res_1);
        EXPECT_EQ(mc_same.size(), 4u);
    }

    // Case 3: residues from different chains
    {
        const gemmi::Chain* chain_B = nullptr;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "B") {
                chain_B = &chain;
                break;
            }
        }
        if (chain_B) {
            if (!chain_B->residues.empty()) {
                const gemmi::Residue& res_b_1 = chain_B->residues[0];
                std::vector<const gemmi::Atom*> mc_diff = coot::crankshaft::get_mainchain_atoms_gemmi(res_1, res_b_1);
                EXPECT_EQ(mc_diff.size(), 4u);
            }
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}