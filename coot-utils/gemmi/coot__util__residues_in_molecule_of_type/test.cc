#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, ResiduesInMoleculeOfType) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // case 1: typical residues (e.g., ALA should exist in the PDB)
    {
        std::string res_type = "ALA";
        std::vector<gemmi::CRA> v = coot::util::residues_in_molecule_of_type_gemmi(st, res_type);
        EXPECT_EQ(v.size(), 18);
        if (v.size() > 0) {
            EXPECT_EQ(v[0].chain->name, "A");
            EXPECT_EQ(v[0].residue->seqid.num.value, 21);
            EXPECT_EQ(v[0].residue->name, "ALA");
        }
    }

    // case 2: non-existent residue type (edge case)
    {
        std::string res_type = "ZZZ";
        std::vector<gemmi::CRA> v = coot::util::residues_in_molecule_of_type_gemmi(st, res_type);
        EXPECT_EQ(v.size(), 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}