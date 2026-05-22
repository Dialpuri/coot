#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_penultimate_residue_in_chain) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // case 1: chain A has residues, get penultimate
    {
        const gemmi::Model& model = st.models[0];
        auto it = std::find_if(model.chains.begin(), model.chains.end(),
                               [](const gemmi::Chain& c) { return c.name == "A"; });
        EXPECT_NE(it, model.chains.end()) << "chain should exist";
        
        size_t nres = it->residues.size();
        EXPECT_EQ(nres, 267);
        
        const gemmi::Residue* penultimate = coot::util::get_penultimate_residue_in_chain_gemmi(*it);
        EXPECT_NE(penultimate, nullptr) << "penultimate residue should be found";
        if (penultimate) {
            EXPECT_EQ(penultimate->seqid.num.value, 297);
            EXPECT_EQ(penultimate->name, "ARG");
        }
    }
    
    // case 2: chain Z does not exist
    {
        const gemmi::Model& model = st.models[0];
        auto it = std::find_if(model.chains.begin(), model.chains.end(),
                               [](const gemmi::Chain& c) { return c.name == "Z"; });
        EXPECT_EQ(it, model.chains.end()) << "chain should not exist";
        
        // Pass a null chain (we can't construct one, so just don't call the function)
        // or test with an empty chain
        const gemmi::Chain empty_chain;
        const gemmi::Residue* penultimate = coot::util::get_penultimate_residue_in_chain_gemmi(empty_chain);
        EXPECT_EQ(penultimate, nullptr) << "penultimate should be nullptr for empty chain";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}