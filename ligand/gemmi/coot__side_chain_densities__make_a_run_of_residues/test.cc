#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, SideChainDensitiesMakeARunOfResidues) {
    // Load PDB using gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Case 1: Valid chain and residue range that should exist
    {
        coot::side_chain_densities scd;
        std::string chain_id = "A";
        int resno_start = 10;
        int resno_end = 20;
        
        std::vector<gemmi::CRA> residues = scd.make_a_run_of_residues_gemmi(st, chain_id, resno_start, resno_end);
        
        EXPECT_EQ(residues.size(), 11u);
        
        // Print individual residues
        for (size_t i = 0; i < residues.size(); i++) {
            gemmi::Residue* res = residues[i].residue;
            int seq_num = res->seqid.num.value;
            std::string res_name = res->name;
            // Expected values from oracle
            if (i == 0) { EXPECT_EQ(seq_num, 10); EXPECT_EQ(res_name, "ILE"); }
            else if (i == 1) { EXPECT_EQ(seq_num, 11); EXPECT_EQ(res_name, "GLY"); }
            else if (i == 2) { EXPECT_EQ(seq_num, 12); EXPECT_EQ(res_name, "GLU"); }
            else if (i == 3) { EXPECT_EQ(seq_num, 13); EXPECT_EQ(res_name, "GLY"); }
            else if (i == 4) { EXPECT_EQ(seq_num, 14); EXPECT_EQ(res_name, "THR"); }
            else if (i == 5) { EXPECT_EQ(seq_num, 15); EXPECT_EQ(res_name, "TYR"); }
            else if (i == 6) { EXPECT_EQ(seq_num, 16); EXPECT_EQ(res_name, "GLY"); }
            else if (i == 7) { EXPECT_EQ(seq_num, 17); EXPECT_EQ(res_name, "VAL"); }
            else if (i == 8) { EXPECT_EQ(seq_num, 18); EXPECT_EQ(res_name, "VAL"); }
            else if (i == 9) { EXPECT_EQ(seq_num, 19); EXPECT_EQ(res_name, "TYR"); }
            else if (i == 10) { EXPECT_EQ(seq_num, 20); EXPECT_EQ(res_name, "LYS"); }
        }
    }
    
    // Case 2: Non-existent chain
    {
        coot::side_chain_densities scd;
        std::string chain_id = "Z";
        int resno_start = 10;
        int resno_end = 20;
        
        std::vector<gemmi::CRA> residues = scd.make_a_run_of_residues_gemmi(st, chain_id, resno_start, resno_end);
        
        EXPECT_EQ(residues.size(), 0u);
    }
    
    // Case 3: Valid chain but non-existent residue range
    {
        coot::side_chain_densities scd;
        std::string chain_id = "A";
        int resno_start = 999;
        int resno_end = 1000;
        
        std::vector<gemmi::CRA> residues = scd.make_a_run_of_residues_gemmi(st, chain_id, resno_start, resno_end);
        
        EXPECT_EQ(residues.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}