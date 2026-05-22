#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper function to count atoms in a structure
static int count_atoms(const gemmi::Structure& st) {
    int atom_count = 0;
    for (const auto& model : st.models) {
        for (const auto& chain : model.chains) {
            for (const auto& residue : chain.residues) {
                atom_count += residue.atoms.size();
            }
        }
    }
    return atom_count;
}

TEST(OracleTest, TrimToMainchain) {
    // Test case 1: First call to trim_to_mainchain
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
        int initial_atoms = count_atoms(st);
        EXPECT_EQ(initial_atoms, 4280);

        coot::db_strands::trim_to_mainchain_gemmi(st);

        int final_atoms = count_atoms(st);
        EXPECT_EQ(final_atoms, 1323);
        EXPECT_EQ(initial_atoms - final_atoms, 2957);
    }

    // Test case 2: Second call on already trimmed molecule (should remove 0 atoms)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
        
        coot::db_strands::trim_to_mainchain_gemmi(st);
        
        int before_second = count_atoms(st);
        EXPECT_EQ(before_second, 1323);

        coot::db_strands::trim_to_mainchain_gemmi(st);
        
        int after_second = count_atoms(st);
        EXPECT_EQ(after_second, 1323);
        EXPECT_EQ(before_second - after_second, 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}