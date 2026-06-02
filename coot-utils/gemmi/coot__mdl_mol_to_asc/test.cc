#include <gtest/gtest.h>
#include <fstream>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, mdl_mol_to_asc_gemmi) {
    // Create a minimal V2000 mol file (water molecule)
    std::string mol_content =
        "Water\n"
        "  RDKit\n"
        "  3D\n"
        "\n"
        "  3  2  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    0.7570    0.5866    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "   -0.7570    0.5866    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  1  2  1  0  0  0  0\n"
        "  1  3  1  0  0  0  0\n"
        "M  END\n";

    // Write to a temp file
    std::string mol_path = "/tmp/test_water.mol";
    {
        std::ofstream ofs(mol_path);
        ofs << mol_content;
        ofs.close();
    }

    // Case 1: valid mol file
    {
        lig_build::molfile_molecule_t m;
        m.read(mol_path);

        EXPECT_EQ(m.atoms.size(), 0u);  // oracle: mol_atoms: 0

        auto asc = coot::mdl_mol_to_asc_gemmi(m);

        EXPECT_EQ(asc.read_success, 0);   // oracle: read_success: 0
        EXPECT_EQ(asc.n_selected_atoms, 0); // oracle: n_selected_atoms: 0
        EXPECT_EQ(asc.mol, nullptr);       // oracle: mol_nullptr: true

        if (asc.mol) delete asc.mol;
    }

    // Case 2: empty mol (no file read)
    {
        lig_build::molfile_molecule_t m;
        // Don't read anything — should produce empty/failed result

        EXPECT_EQ(m.atoms.size(), 0u);  // oracle: mol_atoms: 0

        auto asc = coot::mdl_mol_to_asc_gemmi(m);

        EXPECT_EQ(asc.read_success, 0);   // oracle: read_success: 0
        EXPECT_EQ(asc.n_selected_atoms, 0); // oracle: n_selected_atoms: 0
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}