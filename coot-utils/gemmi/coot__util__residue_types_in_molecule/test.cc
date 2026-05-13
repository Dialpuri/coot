#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, ResidueTypesInMolecule) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    {
        std::vector<std::string> res_types = coot::util::residue_types_in_molecule_gemmi(st);
        std::string res_types_str;
        for (const auto &s : res_types) res_types_str += " " + s;
        
        EXPECT_EQ(res_types.size(), 21);
        EXPECT_EQ(res_types_str, " MET GLU ASN PHE GLN LYS VAL ILE GLY THR TYR ALA ARG LEU SER HIS PRO ASP CYS TRP HOH");
    }

    {
        gemmi::Structure empty_st;
        std::vector<std::string> res_types = coot::util::residue_types_in_molecule_gemmi(empty_st);
        EXPECT_EQ(res_types.size(), 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}