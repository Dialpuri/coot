#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, restraints_container_t_get_atom_index) {
    // Load the PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Get a real residue: chain A, residue 1 (MET)
    const gemmi::Model &model = st.models[0];

    const gemmi::Chain *chain_p = nullptr;
    for (const auto &chain : model.chains) {
        if (chain.name == "A") {
            chain_p = &chain;
            break;
        }
    }
    ASSERT_NE(chain_p, nullptr);

    const gemmi::Residue *residue_p = nullptr;
    for (const auto &res : chain_p->residues) {
        if (res.seqid.num.value == 1) {
            residue_p = &res;
            break;
        }
    }
    ASSERT_NE(residue_p, nullptr);

    // Case 1: valid atom name " N  " (N in residue A/1 MET)
    {
        std::string atom_name_in(" N  ");
        int idx = coot::get_atom_index_gemmi(atom_name_in, *residue_p);
        EXPECT_EQ(idx, 0);
    }

    // Case 2: valid atom name " CA " (CA in residue A/1 MET)
    {
        std::string atom_name_in(" CA ");
        int idx = coot::get_atom_index_gemmi(atom_name_in, *residue_p);
        EXPECT_EQ(idx, 1);
    }

    // Case 3: invalid atom name — should return -2
    {
        std::string atom_name_in(" ZZ ");
        int idx = coot::get_atom_index_gemmi(atom_name_in, *residue_p);
        EXPECT_EQ(idx, -2);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}