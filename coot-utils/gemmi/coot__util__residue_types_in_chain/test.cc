#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, ResidueTypesInChain) {
    // case 1: valid chain — get residue types from chain A
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        
        const gemmi::Model &model = st.models[0];
        const gemmi::Chain &chain = model.chains[0];
        
        std::vector<std::string> v = coot::util::residue_types_in_chain_gemmi(chain);
        
        EXPECT_EQ(chain.name, "A");
        EXPECT_EQ(v.size(), 20U);
        
        // Check specific residue types (alphabetically sorted)
        EXPECT_EQ(v[0], "ALA");
        EXPECT_EQ(v[1], "ARG");
        EXPECT_EQ(v[2], "ASN");
        EXPECT_EQ(v[3], "ASP");
        EXPECT_EQ(v[4], "CYS");
        EXPECT_EQ(v[5], "GLN");
        EXPECT_EQ(v[6], "GLU");
        EXPECT_EQ(v[7], "GLY");
        EXPECT_EQ(v[8], "HIS");
        EXPECT_EQ(v[9], "ILE");
        EXPECT_EQ(v[10], "LEU");
        EXPECT_EQ(v[11], "LYS");
        EXPECT_EQ(v[12], "MET");
        EXPECT_EQ(v[13], "PHE");
        EXPECT_EQ(v[14], "PRO");
        EXPECT_EQ(v[15], "SER");
        EXPECT_EQ(v[16], "THR");
        EXPECT_EQ(v[17], "TRP");
        EXPECT_EQ(v[18], "TYR");
        EXPECT_EQ(v[19], "VAL");
    }

    // case 2: solvent chain — should return empty
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        
        const gemmi::Model &model = st.models[0];
        const gemmi::Chain &chain = model.chains.back(); // last chain
        
        std::vector<std::string> v = coot::util::residue_types_in_chain_gemmi(chain);
        
        EXPECT_EQ(chain.name, "B");
        EXPECT_EQ(v.size(), 1U);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}