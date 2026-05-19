#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, atom_spec_t_get_atom_gemmi) {
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Case 1: valid atom — chain A, residue 10, atom N
    {
        coot::atom_spec_t spec("A", 10, "", " N  ", "");
        const gemmi::Atom *atom = spec.get_atom_gemmi(st);
        EXPECT_TRUE(atom != nullptr);
        EXPECT_EQ(atom->name, "N");
    }

    // Case 2: invalid chain ID — should return null
    {
        coot::atom_spec_t spec("Z", 10, "", " N  ", "");
        const gemmi::Atom *atom = spec.get_atom_gemmi(st);
        EXPECT_EQ(atom, nullptr);
    }

    // Case 3: non-existent residue number — should return null
    {
        coot::atom_spec_t spec("A", 9999, "", " N  ", "");
        const gemmi::Atom *atom = spec.get_atom_gemmi(st);
        EXPECT_EQ(atom, nullptr);
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}