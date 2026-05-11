#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, deep_copy_residue) {
    // Case 1: valid residue - get residue 10 from chain A
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        gemmi::Model& model = st.models[0];
        
        gemmi::Residue* original_res = nullptr;
        for (auto& chain : model.chains)
            for (auto& res : chain.residues)
                if (res.seqid.num.value == 10 && std::string(1, res.seqid.icode) == " ") {
                    original_res = &res;
                    break;
                }
        
        EXPECT_NE(original_res, nullptr) << "residue not found";
        
        if (original_res) {
            // Build CRA for the original residue
            gemmi::Chain* chain_ptr = nullptr;
            for (auto& chain : model.chains)
                for (auto& res : chain.residues)
                    if (&res == original_res) {
                        chain_ptr = &chain;
                        break;
                    }
            gemmi::CRA cra{chain_ptr, original_res, nullptr};
            
            // Print original residue info
            std::string resname = original_res->name;
            int seqnum = original_res->seqid.num.value;
            char inscode = original_res->seqid.icode;
            std::string chainid = chain_ptr->name;
            
            // Make the copy
            gemmi::Residue res_copy = coot::deep_copy_residue_gemmi(cra);
            
            EXPECT_EQ(res_copy.name, "ILE");
            EXPECT_EQ(res_copy.seqid.num.value, 10);
            EXPECT_EQ(res_copy.seqid.icode, ' ');
            
            // Verify atom count
            EXPECT_EQ(res_copy.atoms.size(), 8);
            EXPECT_EQ(original_res->atoms.size(), 8);
        }
    }
    
    // Case 2: another valid residue - test with residue 20
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        gemmi::Model& model = st.models[0];
        
        gemmi::Residue* original_res = nullptr;
        for (auto& chain : model.chains)
            for (auto& res : chain.residues)
                if (res.seqid.num.value == 20 && std::string(1, res.seqid.icode) == " ") {
                    original_res = &res;
                    break;
                }
        
        EXPECT_NE(original_res, nullptr) << "residue not found";
        
        if (original_res) {
            // Build CRA for the original residue
            gemmi::Chain* chain_ptr = nullptr;
            for (auto& chain : model.chains)
                for (auto& res : chain.residues)
                    if (&res == original_res) {
                        chain_ptr = &chain;
                        break;
                    }
            gemmi::CRA cra{chain_ptr, original_res, nullptr};
            
            // Make the copy
            gemmi::Residue res_copy = coot::deep_copy_residue_gemmi(cra);
            
            EXPECT_EQ(res_copy.name, "LYS");
            EXPECT_EQ(res_copy.seqid.num.value, 20);
            EXPECT_EQ(res_copy.seqid.icode, ' ');
            
            // Verify atom count
            EXPECT_EQ(res_copy.atoms.size(), 5);
            EXPECT_EQ(original_res->atoms.size(), 5);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}