#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Parse a CID string like "//A/1/CA" into chain, residue seqid, and atom name
static gemmi::CRA cid_to_atom_gemmi(gemmi::Model& model, const std::string& cid) {
    // Skip the leading "//"
    std::string rest = cid.substr(2);
    
    // Parse chain
    size_t pos = rest.find('/');
    std::string chain_id = rest.substr(0, pos);
    rest = rest.substr(pos + 1);
    
    // Parse residue sequence number
    pos = rest.find('/');
    int res_seq = std::stoi(rest.substr(0, pos));
    rest = rest.substr(pos + 1);
    
    // Remaining is atom name
    std::string atom_name = rest;
    
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == chain_id) {
            for (gemmi::Residue& residue : chain.residues) {
                if (residue.seqid.num.value == res_seq) {
                    for (gemmi::Atom& atom : residue.atoms) {
                        // Trim whitespace for comparison
                        std::string an = atom.name;
                        an.erase(an.find_last_not_of(" \t\r\n") + 1);
                        an.erase(0, an.find_first_not_of(" \t\r\n"));
                        if (an == atom_name) {
                            return gemmi::CRA{&chain, &residue, &atom};
                        }
                    }
                }
            }
        }
    }
    
    return gemmi::CRA{};
}

TEST(OracleTest, is_main_chain_p) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    
    // case 1: main chain CA atom (should be true)
    {
        std::string cid = "//A/1/CA";
        gemmi::CRA cra = cid_to_atom_gemmi(st.models[0], cid);
        ASSERT_NE(cra.atom, nullptr);
        bool result = coot::is_main_chain_p_gemmi(cra);
        EXPECT_TRUE(result);
    }

    // case 2: side chain CB atom (should be false)
    {
        std::string cid = "//A/1/CB";
        gemmi::CRA cra = cid_to_atom_gemmi(st.models[0], cid);
        ASSERT_NE(cra.atom, nullptr);
        bool result = coot::is_main_chain_p_gemmi(cra);
        EXPECT_FALSE(result);
    }

    // case 3: GLY HA2 (special case, should be true)
    {
        std::string cid = "//A/11/HA2";
        gemmi::CRA cra = cid_to_atom_gemmi(st.models[0], cid);
        ASSERT_NE(cra.atom, nullptr);
        bool result = coot::is_main_chain_p_gemmi(cra);
        EXPECT_TRUE(result);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
