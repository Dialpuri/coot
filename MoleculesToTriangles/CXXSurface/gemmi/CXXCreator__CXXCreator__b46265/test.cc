#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, CXXCreator_Selection) {
    // Case 1: Valid Manager with loaded PDB
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        // Count total atoms in the structure for reference
        int total_atoms = CXXCreator::count_all_atoms(st);
        EXPECT_EQ(total_atoms, 2108);

        // Replicate CXXCreator::CXXCreator(mmdb::Manager*) behavior
        std::vector<gemmi::Atom*> SelAtom;
        int nSelAtoms = 0;
        int selHnd = CXXCreator::CXXCreator_gemmi(st, SelAtom, nSelAtoms);

        EXPECT_EQ(selHnd, 1);
        EXPECT_EQ(nSelAtoms, 2107);
        EXPECT_FALSE(nSelAtoms == total_atoms); // sel_atoms_match_total: false
        if (!SelAtom.empty() && nSelAtoms > 0) {
            EXPECT_EQ(CXXCreator::pad_atom_name(*SelAtom[0]), " N  ");
        } else {
            FAIL() << "SelAtom is empty or nSelAtoms is 0";
        }
    }

    // Case 2: Empty Manager (no PDB loaded) — edge case
    {
        gemmi::Structure empty_st;

        int empty_total = 0;
        if (!empty_st.models.empty()) {
            const gemmi::Model& emodel = empty_st.models[0];
            for (const gemmi::Chain& ec : emodel.chains) {
                for (const gemmi::Residue& er : ec.residues) {
                    empty_total += static_cast<int>(er.atoms.size());
                }
            }
        }
        EXPECT_EQ(empty_total, 0);

        // Same constructor logic on empty manager
        std::vector<gemmi::Atom*> SelAtom;
        int nSelAtoms = 0;
        int selHnd = CXXCreator::CXXCreator_gemmi(empty_st, SelAtom, nSelAtoms);

        EXPECT_EQ(selHnd, 1);
        EXPECT_EQ(nSelAtoms, 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}