#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// ... TEST(...) blocks here ...

TEST(OracleTest, mdl_mol_to_asc_gemmi) {
    // Case 1: manually construct a molfile_molecule_t with atoms
    {
        lig_build::molfile_molecule_t m;

        lig_build::molfile_atom_t a1(0.0f, 0.0f, 0.0f, "O");
        a1.name = "O  ";
        m.atoms.push_back(a1);

        lig_build::molfile_atom_t a2(0.757f, 0.587f, 0.0f, "H");
        a2.name = "H1 ";
        m.atoms.push_back(a2);

        lig_build::molfile_atom_t a3(-0.757f, 0.587f, 0.0f, "H");
        a3.name = "H2 ";
        m.atoms.push_back(a3);

        float b_factor = 20.0f;

        coot::atom_selection_container_t_gemmi asc = coot::mdl_mol_to_asc_gemmi(m, b_factor);

        EXPECT_EQ(asc.n_selected_atoms, 3);
        EXPECT_NE(asc.mol, nullptr);
        if (asc.mol && asc.n_selected_atoms > 0) {
            std::string atom_name(asc.mol->models[0].chains[0].residues[0].atoms[0].name);
            EXPECT_EQ(atom_name, "O  ");
        }
        if (asc.mol) delete asc.mol;
    }

    // Case 2: empty molecule (guards the if(m.atoms.size()) path)
    {
        lig_build::molfile_molecule_t m_empty;
        float b_factor = 30.0f;

        coot::atom_selection_container_t_gemmi asc = coot::mdl_mol_to_asc_gemmi(m_empty, b_factor);

        EXPECT_EQ(asc.n_selected_atoms, 0);
        EXPECT_EQ(asc.mol, nullptr);
        if (asc.mol) delete asc.mol;
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}