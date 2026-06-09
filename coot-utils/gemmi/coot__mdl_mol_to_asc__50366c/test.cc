#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, mdl_mol_to_asc_gemmi) {
    // Case 1: Build a small molecule manually (ethane-like: 2 carbons + 6 hydrogens)
    {
        lig_build::molfile_molecule_t mol;

        // Add 8 atoms: 2 carbons + 6 hydrogens
        mol.atoms.push_back(lig_build::molfile_atom_t(0.0f, 0.0f, 0.0f, "C"));
        mol.atoms.push_back(lig_build::molfile_atom_t(1.54f, 0.0f, 0.0f, "C"));
        mol.atoms.push_back(lig_build::molfile_atom_t(-0.37f, 1.0f, 0.0f, "H"));
        mol.atoms.push_back(lig_build::molfile_atom_t(-0.37f, -0.5f, 0.87f, "H"));
        mol.atoms.push_back(lig_build::molfile_atom_t(-0.37f, -0.5f, -0.87f, "H"));
        mol.atoms.push_back(lig_build::molfile_atom_t(1.91f, 1.0f, 0.0f, "H"));
        mol.atoms.push_back(lig_build::molfile_atom_t(1.91f, -0.5f, 0.87f, "H"));
        mol.atoms.push_back(lig_build::molfile_atom_t(1.91f, -0.5f, -0.87f, "H"));

        // Add bonds using SINGLE_BOND enum
        mol.bonds.push_back(lig_build::molfile_bond_t(0, 1, lig_build::bond_t::SINGLE_BOND));
        mol.bonds.push_back(lig_build::molfile_bond_t(0, 2, lig_build::bond_t::SINGLE_BOND));
        mol.bonds.push_back(lig_build::molfile_bond_t(0, 3, lig_build::bond_t::SINGLE_BOND));
        mol.bonds.push_back(lig_build::molfile_bond_t(0, 4, lig_build::bond_t::SINGLE_BOND));
        mol.bonds.push_back(lig_build::molfile_bond_t(1, 5, lig_build::bond_t::SINGLE_BOND));
        mol.bonds.push_back(lig_build::molfile_bond_t(1, 6, lig_build::bond_t::SINGLE_BOND));
        mol.bonds.push_back(lig_build::molfile_bond_t(1, 7, lig_build::bond_t::SINGLE_BOND));

        EXPECT_EQ(mol.atoms.size(), 8u);
        EXPECT_EQ(mol.bonds.size(), 7u);

        gemmi::Structure st = coot::mdl_mol_to_asc_gemmi(mol);

        // read_success: structure has models (non-empty)
        EXPECT_TRUE(!st.models.empty());

        // n_selected_atoms: count atoms in structure
        int total_atoms = 0;
        for (auto& model : st.models)
            for (auto& chain : model.chains)
                for (auto& res : chain.residues)
                    total_atoms += static_cast<int>(res.atoms.size());
        EXPECT_EQ(total_atoms, 8);

        EXPECT_EQ(st.models.size(), 1u);
        {
            auto& model = st.models[0];
            EXPECT_EQ(model.chains.size(), 1u);
            {
                auto& chain = model.chains[0];
                EXPECT_EQ(chain.residues.size(), 1u);
                {
                    auto& res = chain.residues[0];
                    EXPECT_EQ(res.atoms.size(), 8u);
                    EXPECT_EQ(std::string(res.name), "UNL");
                }
            }
        }
    }

    // Case 2: Empty molecule (edge case)
    {
        lig_build::molfile_molecule_t empty_mol;

        EXPECT_EQ(empty_mol.atoms.size(), 0u);
        EXPECT_EQ(empty_mol.bonds.size(), 0u);

        gemmi::Structure st = coot::mdl_mol_to_asc_gemmi(empty_mol);

        // read_success: false means empty structure (no models)
        EXPECT_TRUE(st.models.empty());

        // n_selected_atoms: 0 atoms
        int total_atoms = 0;
        for (auto& model : st.models)
            for (auto& chain : model.chains)
                for (auto& res : chain.residues)
                    total_atoms += static_cast<int>(res.atoms.size());
        EXPECT_EQ(total_atoms, 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
