#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, mol_to_asc_rdkit) {
    // Create a minimal valid mol2 file (simple benzene-like ring)
    std::string mol2_content =
"@<TRIPOS>MOLECULE\n"
"BENZENE\n"
" 12  6  0  0  0\n"
"SMALL\n"
"NO_CHARGE\n"
"\n"
"@<TRIPOS>ATOM\n"
"      1    C1      1.2124  -0.7000   0.0000  C.ar      1    UNL     0.0000\n"
"      2    C2      2.4247  -0.7000   0.0000  C.ar      1    UNL     0.0000\n"
"      3    C3      3.0247  -1.8830   0.0000  C.ar      1    UNL     0.0000\n"
"      4    C4      2.4247  -3.0660   0.0000  C.ar      1    UNL     0.0000\n"
"      5    C5      1.2124  -3.0660   0.0000  C.ar      1    UNL     0.0000\n"
"      6    C6      0.6124  -1.8830   0.0000  C.ar      1    UNL     0.0000\n"
"      7    H1      1.5760   0.2000   0.0000  H         1    UNL     0.0000\n"
"      8    H2      2.8012   0.2000   0.0000  H         1    UNL     0.0000\n"
"      9    H3      3.9689  -1.9000   0.0000  H         1    UNL     0.0000\n"
"     10    H4      2.8012  -3.9660   0.0000  H         1    UNL     0.0000\n"
"     11    H5      0.8359  -3.9660   0.0000  H         1    UNL     0.0000\n"
"     12    H6      0.2359  -1.9000   0.0000  H         1    UNL     0.0000\n"
"\n"
"@<TRIPOS>BOND\n"
"      1      1      2  1.50\n"
"      2      2      3  1.50\n"
"      3      3      4  1.50\n"
"      4      4      5  1.50\n"
"      5      5      6  1.50\n"
"      6      6      1  1.50\n"
"      7      1      7  1\n"
"      8      2      8  1\n"
"      9      3      9  1\n"
"     10      4     10  1\n"
"     11      5     11  1\n"
"     12      6     12  1\n"
"\n"
"@<TRIPOS>SUBSTRUCTURE\n"
"1 BENZENE 12 6 0 0 0 SMALL\n"
;

    // Write the mol2 file
    std::string mol2_path = "/tmp/test_benzene.mol2";
    {
        std::ofstream out(mol2_path);
        out << mol2_content;
    }

    // Case 1: Valid .mol2 file
    {
        atom_selection_container_t asc = coot::mol_to_asc_rdkit_gemmi(mol2_path);
        EXPECT_EQ(asc.read_success, 1);
        EXPECT_EQ(asc.n_selected_atoms, 12);
        EXPECT_NE(asc.mol, nullptr);
        if (asc.mol) {
            EXPECT_EQ(static_cast<int>(asc.mol->models.size()), 1);
        }
    }

    // Case 2: Non-existent file
    {
        std::string bad_file = "/tmp/does_not_exist_xyz.mol2";
        atom_selection_container_t asc2 = coot::mol_to_asc_rdkit_gemmi(bad_file);
        EXPECT_EQ(asc2.read_success, 0);
        EXPECT_EQ(asc2.n_selected_atoms, 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}