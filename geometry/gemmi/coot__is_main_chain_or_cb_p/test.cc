#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, is_main_chain_or_cb_p) {
    using namespace coot;
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];
    
    // case 1: main chain atom (CA) — should be true
    {
        std::string chain_name = "A";
        int seqnum = 10;
        std::string atom_name = "CA";
        
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == chain_name) {
                for (gemmi::Residue& res : chain.residues) {
                    if (res.seqid.num.value == seqnum) {
                        for (gemmi::Atom& atom : res.atoms) {
                            if (atom.name == atom_name) {
                                bool result = is_main_chain_or_cb_p_gemmi(atom);
                                EXPECT_TRUE(result);
                                goto next_test;
                            }
                        }
                    }
                }
            }
        }
        next_test:
        (void)0;
    }

    // case 2: sidechain beta carbon (CB) — should be true (per function spec)
    {
        std::string chain_name = "A";
        int seqnum = 10;
        std::string atom_name = "CB";
        
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == chain_name) {
                for (gemmi::Residue& res : chain.residues) {
                    if (res.seqid.num.value == seqnum) {
                        for (gemmi::Atom& atom : res.atoms) {
                            if (atom.name == atom_name) {
                                bool result = is_main_chain_or_cb_p_gemmi(atom);
                                EXPECT_TRUE(result);
                                goto next_test2;
                            }
                        }
                    }
                }
            }
        }
        next_test2:
        (void)0;
    }

    // case 3: sidechain gamma carbon (CG1) in ILE 10 — should be false
    {
        std::string chain_name = "A";
        int seqnum = 10;
        std::string atom_name = "CG1";
        
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == chain_name) {
                for (gemmi::Residue& res : chain.residues) {
                    if (res.seqid.num.value == seqnum) {
                        for (gemmi::Atom& atom : res.atoms) {
                            if (atom.name == atom_name) {
                                bool result = is_main_chain_or_cb_p_gemmi(atom);
                                EXPECT_FALSE(result);
                                goto next_test3;
                            }
                        }
                    }
                }
            }
        }
        next_test3:
        (void)0;
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}