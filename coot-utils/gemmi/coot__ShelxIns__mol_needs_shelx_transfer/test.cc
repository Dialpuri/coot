#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, mol_needs_shelx_transfer) {
    // case: empty structure (equivalent to null molecule)
    {
        coot::ShelxIns shelx;
        gemmi::Structure st;
        bool result = shelx.mol_needs_shelx_transfer_gemmi(st);
        EXPECT_EQ(result, true);
    }

    // case: structure with 2 chains (example.pdb)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        
        coot::ShelxIns shelx;
        bool result = shelx.mol_needs_shelx_transfer_gemmi(st);
        EXPECT_EQ(result, false);
    }

    // case: structure with 1 chain with empty chain ID
    {
        gemmi::Structure st;
        gemmi::Model model;
        gemmi::Chain chain;
        chain.name = "";
        model.chains.push_back(chain);
        
        gemmi::Residue res;
        res.name = "ALA";
        res.seqid.num.value = 1;
        res.seqid.icode = ' ';
        chain.residues.push_back(res);
        
        st.models.push_back(model);
        
        coot::ShelxIns shelx;
        bool result = shelx.mol_needs_shelx_transfer_gemmi(st);
        EXPECT_EQ(result, true);
    }

    // case: structure with 1 chain with non-empty chain ID
    {
        gemmi::Structure st;
        gemmi::Model model;
        gemmi::Chain chain;
        chain.name = "A";
        model.chains.push_back(chain);
        
        gemmi::Residue res;
        res.name = "ALA";
        res.seqid.num.value = 1;
        res.seqid.icode = ' ';
        chain.residues.push_back(res);
        
        st.models.push_back(model);
        
        coot::ShelxIns shelx;
        bool result = shelx.mol_needs_shelx_transfer_gemmi(st);
        EXPECT_EQ(result, false);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}