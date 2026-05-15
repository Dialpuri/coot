#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, db_strands_get_mol_gemmi) {
    // case 1: valid PDB
    {
        std::string filename = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";
        gemmi::Structure *mol = coot::db_strands::get_mol_gemmi(filename);
        EXPECT_FALSE(mol == nullptr);
        if (mol) {
            int nModels = static_cast<int>(mol->models.size());
            EXPECT_EQ(nModels, 1);
            delete mol;
        }
    }

    // case 2: non-existent file
    {
        std::string filename = "/nonexistent/file.pdb";
        gemmi::Structure *mol = coot::db_strands::get_mol_gemmi(filename);
        EXPECT_TRUE(mol == nullptr);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}