#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, ClashableAltConfs) {
    // Setup: load structure
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty()) << "Failed to load structure";
    
    gemmi::Model& model = st.models[0];
    
    // Get first atom (CA at residue 10 in chain A)
    gemmi::Atom* at_1 = nullptr;
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == "A") {
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == 10) {
                    for (gemmi::Atom& atom : res.atoms) {
                        if (atom.name == "CA") {
                            at_1 = &atom;
                            break;
                        }
                    }
                }
                if (at_1) break;
            }
        }
        if (at_1) break;
    }
    ASSERT_NE(at_1, nullptr) << "First atom (CA at residue 10) not found";
    
    // Get second atom (CA at residue 11 in chain A)
    gemmi::Atom* at_2 = nullptr;
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == "A") {
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == 11) {
                    for (gemmi::Atom& atom : res.atoms) {
                        if (atom.name == "CA") {
                            at_2 = &atom;
                            break;
                        }
                    }
                }
                if (at_2) break;
            }
        }
        if (at_2) break;
    }
    ASSERT_NE(at_2, nullptr) << "Second atom (CA at residue 11) not found";
    
    // Test case 1: Both atoms with altloc (empty/space in gemmi means '\0')
    {
        bool result = coot::clashable_alt_confs_gemmi(*at_1, *at_2);
        EXPECT_TRUE(result) << "clashable_alt_confs should return true when both altlocs are empty/\\0";
    }
    
    // Test case 2: Same atom (edge case)
    {
        bool result = coot::clashable_alt_confs_gemmi(*at_1, *at_1);
        EXPECT_TRUE(result) << "clashable_alt_confs should return true for same atom";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}