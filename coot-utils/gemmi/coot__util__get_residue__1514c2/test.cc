#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, GetResidue) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // case 1: valid residue in chain A
    {
        std::string chain_id = "A";
        int reso = 50;
        std::string ins_code = "";
        gemmi::CRA cra = coot::util::get_residue_gemmi(chain_id, reso, ins_code, st);
        EXPECT_TRUE(cra.chain != nullptr && cra.residue != nullptr) << "Residue should be found";
        if (cra.chain && cra.residue) {
            EXPECT_EQ(cra.residue->seqid.num.value, 50);
            EXPECT_EQ(cra.residue->name, "ARG");
            EXPECT_EQ(cra.chain->name, "A");
        }
    }
    
    // case 2: non-existent residue
    {
        std::string chain_id = "A";
        int reso = 9999;
        std::string ins_code = "";
        gemmi::CRA cra = coot::util::get_residue_gemmi(chain_id, reso, ins_code, st);
        EXPECT_EQ(cra.chain, nullptr) << "Residue should not be found";
    }
    
    // case 3: residue with insertion code (if exists)
    {
        std::string chain_id = "A";
        int reso = 46;  // Looking at the PDB, residue 46 exists
        std::string ins_code = "";
        gemmi::CRA cra = coot::util::get_residue_gemmi(chain_id, reso, ins_code, st);
        EXPECT_TRUE(cra.chain != nullptr && cra.residue != nullptr) << "Residue should be found";
        if (cra.chain && cra.residue) {
            EXPECT_EQ(cra.residue->seqid.num.value, 46);
            EXPECT_EQ(cra.residue->name, "SER");
            EXPECT_EQ(cra.residue->seqid.icode, ' ');
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}