#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, delete_aniso_records_from_atoms) {
    {
        // Case: all atoms have aniso flag set, function should clear all of them
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        // Count total atoms and set aniso on each
        int total_atoms = 0;
        for (gemmi::Model &model : st.models) {
            for (gemmi::Chain &chain : model.chains) {
                for (gemmi::Residue &residue : chain.residues) {
                    for (gemmi::Atom &atom : residue.atoms) {
                        total_atoms++;
                        // Set non-zero aniso values to simulate "aniso flag set"
                        atom.aniso.u11 = 100;
                        atom.aniso.u22 = 100;
                        atom.aniso.u33 = 100;
                        atom.aniso.u12 = 0;
                        atom.aniso.u13 = 0;
                        atom.aniso.u23 = 0;
                    }
                }
            }
        }

        // Count atoms with aniso flag BEFORE
        int aniso_before = 0;
        for (gemmi::Model &model : st.models) {
            for (gemmi::Chain &chain : model.chains) {
                for (gemmi::Residue &residue : chain.residues) {
                    for (gemmi::Atom &atom : residue.atoms) {
                        if (atom.aniso.nonzero()) {
                            aniso_before++;
                        }
                    }
                }
            }
        }

        EXPECT_EQ(total_atoms, 2107);
        EXPECT_EQ(aniso_before, 2107);

        // Call the function under test
        coot::delete_aniso_records_from_atoms_gemmi(st);

        // Count atoms with aniso flag AFTER
        int aniso_after = 0;
        for (gemmi::Model &model : st.models) {
            for (gemmi::Chain &chain : model.chains) {
                for (gemmi::Residue &residue : chain.residues) {
                    for (gemmi::Atom &atom : residue.atoms) {
                        if (atom.aniso.nonzero()) {
                            aniso_after++;
                        }
                    }
                }
            }
        }

        EXPECT_EQ(aniso_after, 0);
    }

    {
        // Case: no aniso flags set — function should not crash and aniso count stays 0
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        // Do NOT set any aniso flags

        EXPECT_NO_THROW(coot::delete_aniso_records_from_atoms_gemmi(st));

        // Count atoms with aniso flag AFTER (should still be 0)
        int aniso_after = 0;
        for (gemmi::Model &model : st.models) {
            for (gemmi::Chain &chain : model.chains) {
                for (gemmi::Residue &residue : chain.residues) {
                    for (gemmi::Atom &atom : residue.atoms) {
                        if (atom.aniso.nonzero()) {
                            aniso_after++;
                        }
                    }
                }
            }
        }

        EXPECT_EQ(aniso_after, 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}