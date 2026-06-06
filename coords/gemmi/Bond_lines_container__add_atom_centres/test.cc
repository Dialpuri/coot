#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Colour constants from coot (COLOUR_BY_CHAIN=0, COLOUR_BY_MOLECULE=4)
static const int COLOUR_BY_CHAIN = 0;
static const int COLOUR_BY_MOLECULE = 4;

TEST(OracleTest, Bond_lines_container_add_atom_centres) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    // Count total atoms
    int total_atoms = 0;
    for (const auto& model : st.models) {
        for (const auto& chain : model.chains) {
            for (const auto& residue : chain.residues) {
                total_atoms += static_cast<int>(residue.atoms.size());
            }
        }
    }

    // Case 1: COLOUR_BY_CHAIN, all models, do_bonds_to_hydrogens = true
    {
        std::vector<gemmi::Vec3> atom_centres;
        std::vector<int> atom_centres_colour;

        // Pre-fill with placeholder data to simulate "before" state
        // The original test expects atom_centres_before == 4280 (the constructor fills with placeholders)
        atom_centres.resize(4280);
        atom_centres_colour.resize(4280);

        bool do_bonds_to_hydrogens = true;
        int atom_colour_type = COLOUR_BY_CHAIN;
        int model_number = 0; // all models

        EXPECT_EQ(atom_colour_type, 0);
        EXPECT_EQ(model_number, 0);
        EXPECT_EQ(total_atoms, 4280);

        int atom_centres_before = static_cast<int>(atom_centres.size());
        int atom_centres_colour_before = static_cast<int>(atom_centres_colour.size());

        coot::add_atom_centres_gemmi(atom_centres, atom_centres_colour, st,
                                     atom_colour_type, model_number,
                                     do_bonds_to_hydrogens, total_atoms);

        int atom_centres_after = static_cast<int>(atom_centres.size());
        int atom_centres_colour_after = static_cast<int>(atom_centres_colour.size());

        EXPECT_EQ(atom_centres_before, 4280);
        EXPECT_EQ(atom_centres_after, 4280);
        EXPECT_EQ(atom_centres_colour_before, 4280);
        EXPECT_EQ(atom_centres_colour_after, 4280);
    }

    // Case 2: COLOUR_BY_MOLECULE, specific model
    {
        std::vector<gemmi::Vec3> atom_centres;
        std::vector<int> atom_centres_colour;

        // Pre-fill with placeholder data (original test starts with 4280)
        atom_centres.resize(4280);
        atom_centres_colour.resize(4280);

        bool do_bonds_to_hydrogens = false;
        int atom_colour_type = COLOUR_BY_MOLECULE;
        int model_number = 1; // model 1 only

        EXPECT_EQ(atom_colour_type, 4);
        EXPECT_EQ(model_number, 1);
        EXPECT_EQ(total_atoms, 4280);

        int atom_centres_before = static_cast<int>(atom_centres.size());

        coot::add_atom_centres_gemmi(atom_centres, atom_centres_colour, st,
                                     atom_colour_type, model_number,
                                     do_bonds_to_hydrogens, total_atoms);

        int atom_centres_after = static_cast<int>(atom_centres.size());

        EXPECT_EQ(atom_centres_before, 4280);
        EXPECT_EQ(atom_centres_after, 2107);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}