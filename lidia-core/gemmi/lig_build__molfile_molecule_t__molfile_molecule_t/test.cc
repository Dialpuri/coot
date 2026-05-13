#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, molfile_molecule_t) {
    // Test 1: Default constructor creates empty molecule
    {
        lig_build::molfile_molecule_t mol;
        EXPECT_EQ(mol.atoms.size(), 0u);
        EXPECT_EQ(mol.bonds.size(), 0u);
    }

    // Test 2: read() on unsupported file format (PDB) — shows guard behavior
    {
        std::string pdb_path = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";
        lig_build::molfile_molecule_t mol;
        mol.read(pdb_path);
        EXPECT_EQ(mol.atoms.size(), 0u);
        EXPECT_EQ(mol.bonds.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}