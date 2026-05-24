#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, CisTransConversion) {
    // Set up the molecules container and load the test PDB
    gemmi::Structure mol = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Structure standard_residues_mol;

    // Get two consecutive residues (e.g., ALA and CYS at positions 9 and 10 in chain A)
    gemmi::Residue *res1 = nullptr;
    gemmi::Residue *res2 = nullptr;
    
    for (auto& model : mol.models) {
        for (auto& chain : model.chains) {
            for (auto& residue : chain.residues) {
                if (residue.seqid.num.value == 9 && res1 == nullptr) {
                    res1 = &residue;
                }
                if (residue.seqid.num.value == 10) {
                    res2 = &residue;
                }
            }
        }
    }
    
    ASSERT_NE(res1, nullptr) << "Could not get residue 1 (position 9)";
    ASSERT_NE(res2, nullptr) << "Could not get residue 2 (position 10)";
    
    // Case: empty standard_residues_mol (no TNS/CIS residues)
    // Based on oracle: INPUT res1: //A/9, INPUT res2: //A/10, OUTPUT conversion_count: 0
    int result_empty = coot::util::cis_trans_conversion_gemmi(res1, res2, mol, standard_residues_mol);
    EXPECT_EQ(result_empty, 0) << "Conversion count should be 0 with empty standard residues molecule";
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}