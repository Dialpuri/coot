#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_types_in_molecule) {
    // INPUT: PDB path
    std::string pdb_path = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";

    // Read the PDB file using gemmi
    gemmi::Structure st = gemmi::read_pdb_file(pdb_path);

    // Call the function under observation
    std::vector<std::string> types = coot::get_types_in_molecule_gemmi(st);

    // OUTPUT: the returned vector (size + elements)
    EXPECT_EQ(types.size(), 21);
    EXPECT_EQ(types[0], "ALA");
    EXPECT_EQ(types[1], "ARG");
    EXPECT_EQ(types[2], "ASN");
    EXPECT_EQ(types[3], "ASP");
    EXPECT_EQ(types[4], "CYS");
    EXPECT_EQ(types[5], "GLN");
    EXPECT_EQ(types[6], "GLU");
    EXPECT_EQ(types[7], "GLY");
    EXPECT_EQ(types[8], "HIS");
    EXPECT_EQ(types[9], "HOH");
    EXPECT_EQ(types[10], "ILE");
    EXPECT_EQ(types[11], "LEU");
    EXPECT_EQ(types[12], "LYS");
    EXPECT_EQ(types[13], "MET");
    EXPECT_EQ(types[14], "PHE");
    EXPECT_EQ(types[15], "PRO");
    EXPECT_EQ(types[16], "SER");
    EXPECT_EQ(types[17], "THR");
    EXPECT_EQ(types[18], "TRP");
    EXPECT_EQ(types[19], "TYR");
    EXPECT_EQ(types[20], "VAL");
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}