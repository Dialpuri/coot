#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_nth_residue) {
    // Load the example PDB file
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Case 1: nth = 1 - should return first residue (first residue of chain A)
    {
        int nth = 1;
        gemmi::CRA cra = coot::util::get_nth_residue_gemmi(nth, st);
        EXPECT_TRUE(cra.residue != nullptr);
        if (cra.residue) {
            EXPECT_EQ(cra.chain->name, "A");
            EXPECT_EQ(cra.residue->seqid.num.value, 1);
            EXPECT_EQ(cra.residue->name, "MET");
        }
    }
    
    // Case 2: nth = 2 - should return second residue
    {
        int nth = 2;
        gemmi::CRA cra = coot::util::get_nth_residue_gemmi(nth, st);
        EXPECT_TRUE(cra.residue != nullptr);
        if (cra.residue) {
            EXPECT_EQ(cra.chain->name, "A");
            EXPECT_EQ(cra.residue->seqid.num.value, 2);
            EXPECT_EQ(cra.residue->name, "GLU");
        }
    }
    
    // Case 3: nth = 0 - should return nullptr per doc
    {
        int nth = 0;
        gemmi::CRA cra = coot::util::get_nth_residue_gemmi(nth, st);
        EXPECT_FALSE(cra.residue != nullptr);
    }
    
    // Case 4: nth = 1000 - beyond available residues, should return nullptr
    {
        int nth = 1000;
        gemmi::CRA cra = coot::util::get_nth_residue_gemmi(nth, st);
        EXPECT_FALSE(cra.residue != nullptr);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}