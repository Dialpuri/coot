#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, ResiduesWithInsertionCodes) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    int total_residues = 0;
    for (const gemmi::Chain& chain : st.models[0].chains) {
        total_residues += chain.residues.size();
    }
    
    std::vector<gemmi::Residue*> residues = coot::util::residues_with_insertion_codes_gemmi(st);
    
    // Ground truth from oracle
    EXPECT_EQ(total_residues, 268);
    EXPECT_EQ(residues.size(), 0);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}