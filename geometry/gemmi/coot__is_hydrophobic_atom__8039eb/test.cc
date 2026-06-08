#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, is_hydrophobic_atom) {
    // Case 1: hydrophobic residue (ALA) with non-N/O atom → true
    {
        std::string res_name = "ALA";
        std::string atom_name = " CA";
        bool result = coot::is_hydrophobic_atom_gemmi(res_name, atom_name);
        EXPECT_EQ(result, true);
    }

    // Case 2: hydrophilic residue (ASP) with CB atom → false
    {
        std::string res_name = "ASP";
        std::string atom_name = " CB";
        bool result = coot::is_hydrophobic_atom_gemmi(res_name, atom_name);
        EXPECT_EQ(result, false);
    }

    // Case 3: hydrophobic residue (VAL) but N-containing atom name → false
    {
        std::string res_name = "VAL";
        std::string atom_name = " ND1";
        bool result = coot::is_hydrophobic_atom_gemmi(res_name, atom_name);
        EXPECT_EQ(result, false);
    }

    // Case 4: hydrophobic residue (PHE) with O-containing atom name → false
    {
        std::string res_name = "PHE";
        std::string atom_name = " OE1";
        bool result = coot::is_hydrophobic_atom_gemmi(res_name, atom_name);
        EXPECT_EQ(result, false);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}