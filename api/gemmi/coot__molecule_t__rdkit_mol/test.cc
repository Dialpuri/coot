#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, rdkit_mol) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    gemmi::Model& model = st.models[0];

    // Case 1: valid CID (//A/1 - GLY residue in chain A)
    {
        std::string cid = "//A/1";
        void* mol = coot::molecule_t::rdkit_mol_gemmi(cid, model);
        EXPECT_EQ(mol, nullptr);
    }

    // Case 2: invalid CID
    {
        std::string cid = "//B/9999";
        void* mol = coot::molecule_t::rdkit_mol_gemmi(cid, model);
        EXPECT_EQ(mol, nullptr);
    }

    // Case 3: empty CID
    {
        std::string cid = "";
        void* mol = coot::molecule_t::rdkit_mol_gemmi(cid, model);
        EXPECT_EQ(mol, nullptr);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}