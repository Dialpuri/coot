#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Test that creating a substructure from atom selection works correctly
TEST(OracleTest, CreateMMDBManagerFromAtomSelectionStraight) {
    // case: all atoms selected
    {
        gemmi::Structure orig_st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        
        // Count total atoms in original structure
        int n_atoms_total = 0;
        for (const gemmi::Model& model : orig_st.models) {
            for (const gemmi::Chain& chain : model.chains) {
                for (const gemmi::Residue& res : chain.residues) {
                    for (const gemmi::Atom& atom : res.atoms) {
                        n_atoms_total++;
                    }
                }
            }
        }
        
        // Select all atoms
        auto all_atoms_selector = [](const gemmi::CRA& cra) -> bool {
            return true;
        };
        
        gemmi::Structure sub_st = coot::util::create_mmdbmanager_from_atom_selection_straight_gemmi(
            orig_st, all_atoms_selector);
        
        // Count atoms in substructure
        int n_atoms_sub = 0;
        for (const gemmi::Model& model : sub_st.models) {
            for (const gemmi::Chain& chain : model.chains) {
                for (const gemmi::Residue& res : chain.residues) {
                    for (const gemmi::Atom& atom : res.atoms) {
                        n_atoms_sub++;
                    }
                }
            }
        }
        
        EXPECT_EQ(n_atoms_sub, n_atoms_total) << "substructure atoms";
        EXPECT_EQ(n_atoms_sub, 2107) << "substructure atom count";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}