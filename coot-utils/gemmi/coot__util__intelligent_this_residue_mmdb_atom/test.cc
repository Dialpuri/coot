#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, intelligent_this_residue_mmdb_atom) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // case 1: valid residue with CA atom — expects CA returned
    {
        const gemmi::Chain* chain_a = nullptr;
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                if (chain.name == "A") {
                    chain_a = &chain;
                    break;
                }
            }
            if (chain_a) break;
        }
        
        EXPECT_NE(chain_a, nullptr);
        if (chain_a) {
            const gemmi::Residue* res = nullptr;
            for (const auto& r : chain_a->residues) {
                if (r.seqid.num.value == 10) {
                    res = &r;
                    break;
                }
            }
            EXPECT_NE(res, nullptr);
            if (res) {
                const gemmi::Atom* atom = coot::util::intelligent_this_residue_mmdb_atom_gemmi(*res);
                EXPECT_NE(atom, nullptr);
                if (atom) {
                    EXPECT_TRUE(atom->name.find("CA") != std::string::npos);
                }
            }
        }
    }

    // case 2: valid residue (CA present)
    {
        const gemmi::Chain* chain_a = nullptr;
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                if (chain.name == "A") {
                    chain_a = &chain;
                    break;
                }
            }
            if (chain_a) break;
        }
        
        EXPECT_NE(chain_a, nullptr);
        if (chain_a) {
            const gemmi::Residue* res = nullptr;
            for (const auto& r : chain_a->residues) {
                if (r.seqid.num.value == 1) {
                    res = &r;
                    break;
                }
            }
            EXPECT_NE(res, nullptr);
            if (res) {
                const gemmi::Atom* atom = coot::util::intelligent_this_residue_mmdb_atom_gemmi(*res);
                EXPECT_NE(atom, nullptr);
                if (atom) {
                    EXPECT_TRUE(atom->name.find("CA") != std::string::npos);
                }
            }
        }
    }

    // case 3: invalid chain/residue — should return nullptr
    {
        const gemmi::Chain* chain_x = nullptr;
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                if (chain.name == "X") {
                    chain_x = &chain;
                    break;
                }
            }
            if (chain_x) break;
        }
        EXPECT_EQ(chain_x, nullptr);
        
        const gemmi::Residue* res_9999 = nullptr;
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                for (const auto& r : chain.residues) {
                    if (r.seqid.num.value == 9999) {
                        res_9999 = &r;
                        break;
                    }
                }
            }
            if (res_9999) break;
        }
        EXPECT_EQ(res_9999, nullptr);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}