#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(MinimolMoleculeConstructorGemmi, basic) {
    // Case 1: Single atom with explicit element
    {
        std::vector<gemmi::Vec3> atom_list;
        atom_list.push_back(gemmi::Vec3(1.0, 2.0, 3.0));

        std::string residue_type = "HOH";
        std::string atom_name = "O1";
        std::string chain_id = "X";
        std::string ele = "O";

        coot::minimol::molecule mol = coot::minimol::molecule_gemmi(atom_list, residue_type, atom_name, chain_id, ele);

        EXPECT_EQ(static_cast<int>(mol.fragments.size()), 1);
        EXPECT_EQ(mol.count_atoms(), 1);

        const coot::minimol::fragment &frag = mol.fragments[0];
        EXPECT_EQ(frag.fragment_id, "X");
        EXPECT_EQ(frag[1].name, "HOH");
        EXPECT_EQ(static_cast<int>(frag[1].atoms.size()), 1);

        const coot::minimol::atom &a = frag[1].atoms[0];
        EXPECT_EQ(a.name, "O1");
        EXPECT_EQ(a.element, "O");
        EXPECT_NEAR(a.pos.x(), 1.0, 1e-4);
        EXPECT_NEAR(a.pos.y(), 2.0, 1e-4);
        EXPECT_NEAR(a.pos.z(), 3.0, 1e-4);
    }

    // Case 2: Multiple atoms
    {
        std::vector<gemmi::Vec3> atom_list;
        atom_list.push_back(gemmi::Vec3(0.0, 0.0, 0.0));
        atom_list.push_back(gemmi::Vec3(1.5, 0.0, 0.0));
        atom_list.push_back(gemmi::Vec3(0.75, 1.3, 0.0));

        std::string residue_type = "DUM";
        std::string atom_name = "CA";
        std::string chain_id = "A";
        std::string ele = "C";

        coot::minimol::molecule mol = coot::minimol::molecule_gemmi(atom_list, residue_type, atom_name, chain_id, ele);

        EXPECT_EQ(static_cast<int>(mol.fragments.size()), 1);
        EXPECT_EQ(mol.count_atoms(), 3);
    }

    // Case 3: Empty atom list (edge case)
    {
        std::vector<gemmi::Vec3> atom_list;

        std::string residue_type = "WAT";
        std::string atom_name = "O";
        std::string chain_id = "B";
        std::string ele = "O";

        coot::minimol::molecule mol = coot::minimol::molecule_gemmi(atom_list, residue_type, atom_name, chain_id, ele);

        EXPECT_EQ(static_cast<int>(mol.fragments.size()), 1);
        EXPECT_EQ(mol.count_atoms(), 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
