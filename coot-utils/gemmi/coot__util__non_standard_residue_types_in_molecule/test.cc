#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, NonStandardResidueTypesInMolecule) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // case 1: normal molecule with non-standard residues
    {
        std::vector<std::string> non_standards = coot::util::non_standard_residue_types_in_molecule_gemmi(st);
        EXPECT_EQ(non_standards.size(), 1u);
        EXPECT_EQ(non_standards[0], "HOH");
    }

    // case 2: empty structure (equivalent to null mol guard clause)
    {
        gemmi::Structure empty_st;
        std::vector<std::string> non_standards = coot::util::non_standard_residue_types_in_molecule_gemmi(empty_st);
        EXPECT_EQ(non_standards.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}