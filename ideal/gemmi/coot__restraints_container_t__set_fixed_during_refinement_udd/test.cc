#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// ... TEST(...) blocks here ...

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

TEST(OracleTest, set_fixed_during_refinement_udd) {
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        EXPECT_EQ(st.models.size(), 1);

        gemmi::Model& model = st.models[0];

        // Collect all non-hydrogen atoms from the model
        std::vector<gemmi::Atom*> all_atoms;
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& residue : chain.residues) {
                for (gemmi::Atom& atom : residue.atoms) {
                    if (!atom.element.is_hydrogen()) {
                        all_atoms.push_back(&atom);
                    }
                }
            }
        }

        // Use first 94 atoms (matching the original test's n_atoms)
        constexpr int n_atoms = 94;
        std::vector<gemmi::Atom*> atoms(all_atoms.begin(), all_atoms.begin() + n_atoms);
        EXPECT_EQ(static_cast<int>(atoms.size()), n_atoms);

        // Build fixed_atom_indices: first 82 atoms are fixed (matching original test)
        std::set<int> fixed_atom_indices;
        for (int i = 0; i < 82; i++) {
            fixed_atom_indices.insert(i);
        }

        EXPECT_EQ(static_cast<int>(fixed_atom_indices.size()), 82);

        // Call the target function
        std::vector<bool> fixed_during_refinement;
        coot::set_fixed_during_refinement_udd_gemmi(atoms, n_atoms, fixed_atom_indices, fixed_during_refinement);

        // Verify UDD was set by reading back
        int fixed_count = 0;
        int unfixed_count = 0;

        for (int i = 0; i < n_atoms; i++) {
            if (fixed_during_refinement[i])
                fixed_count++;
            else
                unfixed_count++;
        }

        EXPECT_EQ(fixed_count, 82);
        EXPECT_EQ(unfixed_count, 12);
        EXPECT_EQ(n_atoms, 94);
    }
}