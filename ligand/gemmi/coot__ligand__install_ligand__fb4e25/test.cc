#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Derived class to expose protected members for testing
class testable_ligand : public coot::ligand {
public:
    using coot::ligand::initial_ligand;
};

TEST(OracleTest, install_ligand) {
    testable_ligand lig;

    // Case 1: install first ligand
    {
        EXPECT_EQ(lig.initial_ligand.size(), size_t(0));
        std::string pdb_filename = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";
        coot::install_ligand_gemmi(lig, pdb_filename);
        EXPECT_EQ(lig.initial_ligand.size(), size_t(1));
    }

    // Case 2: install second ligand
    {
        EXPECT_EQ(lig.initial_ligand.size(), size_t(1));
        std::string pdb_filename = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";
        coot::install_ligand_gemmi(lig, pdb_filename);
        EXPECT_EQ(lig.initial_ligand.size(), size_t(2));
    }

    // Case 3: install third ligand
    {
        EXPECT_EQ(lig.initial_ligand.size(), size_t(2));
        std::string pdb_filename = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";
        coot::install_ligand_gemmi(lig, pdb_filename);
        EXPECT_EQ(lig.initial_ligand.size(), size_t(3));
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
