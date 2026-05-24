#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Test for get_simple_molecule
TEST(OracleTest, get_simple_molecule) {
    // Load structure
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    gemmi::Model& model = st.models[0];
    
    // Find residue A/10
    gemmi::Residue* res = nullptr;
    for (auto& chain : model.chains) {
        if (chain.name == "A") {
            for (auto& r : chain.residues) {
                if (r.seqid.num.value == 10) {
                    res = &r;
                    break;
                }
            }
        }
    }
    
    // Case 1: valid residue CID
    {
        bool draw_hydrogens = true;
        coot::simple::molecule_t sm = coot::molecule_t_get_simple_molecule_gemmi(0, res, draw_hydrogens, nullptr);
        EXPECT_TRUE(sm.is_valid());
        EXPECT_EQ(19u, sm.atoms.size());
    }
    
    // Case 2: invalid residue (nullptr)
    {
        bool draw_hydrogens = false;
        coot::simple::molecule_t sm = coot::molecule_t_get_simple_molecule_gemmi(0, nullptr, draw_hydrogens, nullptr);
        EXPECT_FALSE(sm.is_valid());
        EXPECT_EQ(0u, sm.atoms.size());
        EXPECT_EQ(0u, sm.bonds.size());
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}