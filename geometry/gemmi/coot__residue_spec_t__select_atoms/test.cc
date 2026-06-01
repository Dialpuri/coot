#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, residue_spec_t_select_atoms) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Case 1: Valid residue — ARG 50 in chain A (has 11 atoms)
    {
        coot::residue_spec_t spec("A", 50, "");
        auto atoms = coot::select_atoms_gemmi(spec, &st);
        EXPECT_EQ(atoms.size(), 11);
    }

    // Case 2: Non-existent residue — chain A, residue 9999 (should select 0 atoms)
    {
        coot::residue_spec_t spec("A", 9999, "");
        auto atoms = coot::select_atoms_gemmi(spec, &st);
        EXPECT_EQ(atoms.size(), 0);
    }

    // Case 3: null mol — verifies the guard clause (mol == nullptr)
    {
        coot::residue_spec_t spec("A", 1, "");
        auto atoms = coot::select_atoms_gemmi(spec, nullptr);
        EXPECT_EQ(atoms.size(), 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}