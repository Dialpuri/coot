#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, residue_orientation) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty()) << "Failed to load PDB file";

    // Case 1: Valid residue with sidechain atoms (ILE at position 10)
    {
        bool found = false;
        for (gemmi::Model& model : st.models) {
            for (gemmi::Chain& chain : model.chains) {
                for (gemmi::Residue& res : chain.residues) {
                    if (res.seqid.num.value == 10) {
                        found = true;
                        clipper::Mat33<double> identity;
                        clipper::Mat33<double> result = coot::util::residue_orientation_gemmi(res, identity);
                        
                        EXPECT_NEAR(result(0,0), 0.605854, 1e-4);
                        EXPECT_NEAR(result(0,1), 0.174217, 1e-4);
                        EXPECT_NEAR(result(0,2), 0.776267, 1e-4);
                        EXPECT_NEAR(result(1,0), 0.78885, 1e-4);
                        EXPECT_NEAR(result(1,1), -0.0049362, 1e-4);
                        EXPECT_NEAR(result(1,2), -0.614566, 1e-4);
                        EXPECT_NEAR(result(2,0), 0.103236, 1e-4);
                        EXPECT_NEAR(result(2,1), -0.984695, 1e-4);
                        EXPECT_NEAR(result(2,2), 0.140421, 1e-4);
                        EXPECT_NEAR(result.det(), -1.0, 1e-4);
                    }
                }
            }
        }
        EXPECT_TRUE(found) << "Residue ILE at position 10 not found";
    }

    // Case 2: Another valid residue (LYS at position 20)
    {
        bool found = false;
        for (gemmi::Model& model : st.models) {
            for (gemmi::Chain& chain : model.chains) {
                for (gemmi::Residue& res : chain.residues) {
                    if (res.seqid.num.value == 20) {
                        found = true;
                        clipper::Mat33<double> identity;
                        clipper::Mat33<double> result = coot::util::residue_orientation_gemmi(res, identity);
                        
                        EXPECT_NEAR(result(0,0), 0.351098, 1e-4);
                        EXPECT_NEAR(result(0,1), -0.411479, 1e-4);
                        EXPECT_NEAR(result(0,2), -0.84108, 1e-4);
                        EXPECT_NEAR(result(1,0), 0.122246, 1e-4);
                        EXPECT_NEAR(result(1,1), -0.870432, 1e-4);
                        EXPECT_NEAR(result(1,2), 0.476869, 1e-4);
                        EXPECT_NEAR(result(2,0), 0.928324, 1e-4);
                        EXPECT_NEAR(result(2,1), 0.270247, 1e-4);
                        EXPECT_NEAR(result(2,2), 0.255305, 1e-4);
                        EXPECT_NEAR(result.det(), -1.0, 1e-4);
                    }
                }
            }
        }
        EXPECT_TRUE(found) << "Residue LYS at position 20 not found";
    }

    // Case 3: Invalid residue (nullptr case - use empty residue)
    {
        clipper::Mat33<double> identity;
        gemmi::Residue empty_res;
        clipper::Mat33<double> result = coot::util::residue_orientation_gemmi(empty_res, identity);
        
        EXPECT_EQ(result.det(), 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}