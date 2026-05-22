#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// ... TEST(...) blocks here ...

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

// Helper to get test PDB path
static std::string get_test_pdb_path() {
    return "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";
}

TEST(OracleTest, GetHetgroups) {
    // Load the protein-ligand complex which should have hetgroups
    gemmi::Structure st = gemmi::read_pdb_file(get_test_pdb_path());
    
    // Case 1: include waters
    {
        std::vector<std::pair<gemmi::Chain*, gemmi::Residue*>> hetgroups = 
            coot::util::get_hetgroups_gemmi(st, true);
        EXPECT_EQ(hetgroups.size(), 1);
        if (hetgroups.size() > 0) {
            std::pair<gemmi::Chain*, gemmi::Residue*> res_pair = hetgroups[0];
            std::string chain_id = res_pair.first->name;
            std::string res_name = res_pair.second->name;
            EXPECT_EQ(chain_id, "B");
            EXPECT_EQ(res_pair.second->seqid.num.value, 1);
            EXPECT_EQ(res_name, "HOH");
        }
    }
    
    // Case 2: exclude waters
    {
        std::vector<std::pair<gemmi::Chain*, gemmi::Residue*>> hetgroups = 
            coot::util::get_hetgroups_gemmi(st, false);
        EXPECT_EQ(hetgroups.size(), 0);
    }
    
    // Case 3: empty structure (edge case)
    {
        gemmi::Structure empty_st;
        std::vector<std::pair<gemmi::Chain*, gemmi::Residue*>> hetgroups = 
            coot::util::get_hetgroups_gemmi(empty_st, true);
        EXPECT_EQ(hetgroups.size(), 0);
    }
}