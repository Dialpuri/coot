#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, install_ligand) {
    // Case 1: install_ligand with a real gemmi::Structure
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        int n_models = static_cast<int>(st.models.size());
        EXPECT_EQ(n_models, 1);

        coot::ligand lig;

        int size_before = static_cast<int>(lig.initial_ligand.size());
        EXPECT_EQ(size_before, 0);

        lig.install_ligand_gemmi(st);

        int size_after = static_cast<int>(lig.initial_ligand.size());
        EXPECT_EQ(size_after, 1);
    }

    // Case 2: install_ligand twice — verify it keeps growing
    {
        gemmi::Structure st2 = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        coot::ligand lig2;

        int size_before = static_cast<int>(lig2.initial_ligand.size());
        EXPECT_EQ(size_before, 0);

        lig2.install_ligand_gemmi(st2);
        int size_after_first = static_cast<int>(lig2.initial_ligand.size());
        EXPECT_EQ(size_after_first, 1);

        lig2.install_ligand_gemmi(st2);
        int size_after_second = static_cast<int>(lig2.initial_ligand.size());
        EXPECT_EQ(size_after_second, 2);
    }

    // Case 3: install_ligand with a minimal (empty) gemmi::Structure
    {
        gemmi::Structure st3; // empty structure

        coot::ligand lig3;
        int size_before = static_cast<int>(lig3.initial_ligand.size());
        EXPECT_EQ(size_before, 0);

        lig3.install_ligand_gemmi(st3);

        int size_after = static_cast<int>(lig3.initial_ligand.size());
        EXPECT_EQ(size_after, 1);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}