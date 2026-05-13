#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper function to count atoms in a structure
static unsigned int count_atoms(const gemmi::Structure& st) {
    unsigned int count = 0;
    for (const auto& model : st.models) {
        for (const auto& chain : model.chains) {
            for (const auto& residue : chain.residues) {
                count += residue.atoms.size();
            }
        }
    }
    return count;
}

TEST(OracleTest, DeleteLiteralUsingCID) {
    // Load the test PDB file
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    coot::molecule_t mol(std::move(st));
    
    // case 1: delete a single atom (CA of residue 10 in chain A)
    {
        std::string cid = "//A/10/CA";
        unsigned int atoms_before = count_atoms(mol.structure);
        int status = mol.delete_literal_using_cid_gemmi(cid);
        unsigned int atoms_after = count_atoms(mol.structure);

        EXPECT_EQ(status, 1);
        EXPECT_EQ(atoms_before, 2107u);
        EXPECT_EQ(atoms_after, 2106u);
    }

    // case 2: delete a residue range (residues 11 and 12 in chain A)
    {
        std::string cid = "//A/11||//A/12";
        int status = mol.delete_literal_using_cid_gemmi(cid);
        unsigned int atoms_after_range = count_atoms(mol.structure);

        EXPECT_EQ(status, 1);
        EXPECT_EQ(atoms_after_range, 2093u);
    }

    // case 3: invalid cid (should return 0, no atoms deleted)
    {
        std::string cid = "//X/999/CA";
        int status = mol.delete_literal_using_cid_gemmi(cid);
        unsigned int atoms_invalid = count_atoms(mol.structure);

        EXPECT_EQ(status, 0);
        EXPECT_EQ(atoms_invalid, 2093u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}