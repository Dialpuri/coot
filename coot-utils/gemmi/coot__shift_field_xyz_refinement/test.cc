#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, shift_field_xyz_refinement) {
    {
        // Load a PDB file into a gemmi::Structure
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        int n_models = static_cast<int>(st.models.size());
        EXPECT_EQ(n_models, 1);

        // Record starting number of atoms
        int atoms_before = 0;
        if (n_models >= 1) {
            gemmi::Model &model = st.models[0];
            for (gemmi::Chain &chain : model.chains) {
                for (gemmi::Residue &res : chain.residues) {
                    atoms_before += static_cast<int>(res.atoms.size());
                }
            }
        }
        EXPECT_EQ(atoms_before, 2107);

        // Create empty clipper HKL_data objects (no reflections needed since function is a stub)
        clipper::HKL_data<clipper::datatypes::F_sigF<float>> fo0;
        clipper::HKL_data<clipper::data32::Flag> free_flag;

        float resolution = 2.0f;

        // Call the function
        EXPECT_NO_THROW(coot::shift_field_xyz_refinement_gemmi(fo0, free_flag, st, resolution));

        // Record ending number of atoms (should be same since function is a stub)
        int atoms_after = 0;
        if (n_models >= 1) {
            gemmi::Model &model = st.models[0];
            for (gemmi::Chain &chain : model.chains) {
                for (gemmi::Residue &res : chain.residues) {
                    atoms_after += static_cast<int>(res.atoms.size());
                }
            }
        }
        EXPECT_EQ(atoms_after, 2107);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}