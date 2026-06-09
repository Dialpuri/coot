#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(MiniMolMoleculeTest, molecule) {
    // Case 1: Single atom with default element
    {
        std::vector<std::pair<gemmi::Vec3, float>> atom_list;
        atom_list.push_back(std::make_pair(gemmi::Vec3(1.0, 2.0, 3.0), 10.0f));

        std::string residue_type = "ALA";
        std::string atom_name = "CA";
        std::string chain_id = "A";
        std::string ele = "O";

        coot::minimol::molecule mol = coot::minimol::molecule_gemmi(atom_list, residue_type, atom_name, chain_id, ele);

        EXPECT_EQ(mol.fragments.size(), 1u);
        EXPECT_EQ(mol.fragments[0].residues.size(), 2u);
        if (mol.fragments[0].residues.size() > 1) {
            auto &res = mol.fragments[0].residues[1];
            EXPECT_EQ(res.seqnum, 1);
            EXPECT_EQ(res.name, "ALA");
            EXPECT_EQ(res.atoms.size(), 1u);
            if (!res.atoms.empty()) {
                EXPECT_EQ(res.atoms[0].name, "CA");
                EXPECT_EQ(res.atoms[0].element, "O");
                EXPECT_NEAR(res.atoms[0].pos.x(), 1.0, 1e-4);
                EXPECT_NEAR(res.atoms[0].pos.y(), 2.0, 1e-4);
                EXPECT_NEAR(res.atoms[0].pos.z(), 3.0, 1e-4);
                EXPECT_NEAR(res.atoms[0].temperature_factor, 10.0f, 1e-4);
            }
        }
    }

    // Case 2: Multiple atoms with carbon element
    {
        std::vector<std::pair<gemmi::Vec3, float>> atom_list;
        atom_list.push_back(std::make_pair(gemmi::Vec3(1.0, 2.0, 3.0), 15.0f));
        atom_list.push_back(std::make_pair(gemmi::Vec3(4.0, 5.0, 6.0), 20.0f));
        atom_list.push_back(std::make_pair(gemmi::Vec3(7.0, 8.0, 9.0), 25.0f));

        std::string residue_type = "GLY";
        std::string atom_name = "N";
        std::string chain_id = "B";
        std::string ele = "C";

        coot::minimol::molecule mol = coot::minimol::molecule_gemmi(atom_list, residue_type, atom_name, chain_id, ele);

        EXPECT_EQ(mol.fragments.size(), 1u);
        EXPECT_EQ(mol.fragments[0].residues.size(), 4u);
        for (size_t i = 0; i < mol.fragments[0].residues.size(); i++) {
            auto &res = mol.fragments[0].residues[i];
            if (i == 0) {
                EXPECT_EQ(res.seqnum, 1);
                EXPECT_EQ(res.name, "");
                EXPECT_EQ(res.atoms.size(), 0u);
            } else if (i == 1) {
                EXPECT_EQ(res.seqnum, 1);
                EXPECT_EQ(res.name, "GLY");
                EXPECT_EQ(res.atoms.size(), 1u);
            } else if (i == 2) {
                EXPECT_EQ(res.seqnum, 2);
                EXPECT_EQ(res.name, "GLY");
                EXPECT_EQ(res.atoms.size(), 1u);
            } else if (i == 3) {
                EXPECT_EQ(res.seqnum, 3);
                EXPECT_EQ(res.name, "GLY");
                EXPECT_EQ(res.atoms.size(), 1u);
            }
        }
    }

    // Case 3: Empty atom list (edge case)
    {
        std::vector<std::pair<gemmi::Vec3, float>> atom_list;

        std::string residue_type = "WAT";
        std::string atom_name = "O";
        std::string chain_id = "C";
        std::string ele = "O";

        coot::minimol::molecule mol = coot::minimol::molecule_gemmi(atom_list, residue_type, atom_name, chain_id, ele);

        EXPECT_EQ(mol.fragments.size(), 1u);
        EXPECT_EQ(mol.fragments[0].residues.size(), 1u);
        EXPECT_EQ(mol.fragments[0].fragment_id, "C");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
