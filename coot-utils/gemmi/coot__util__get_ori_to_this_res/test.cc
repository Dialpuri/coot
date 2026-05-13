#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, GetOriToThisRes) {
    // Setup: load the structure
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());

    // Case 1: Valid residue with mainchain atoms (ILE at A/10)
    {
        const gemmi::Model& model = st.models[0];
        const gemmi::Chain* chain_a = nullptr;
        for (const gemmi::Chain& chain : model.chains) {
            if (chain.name == "A") {
                chain_a = &chain;
                break;
            }
        }
        ASSERT_NE(chain_a, nullptr);
        
        const gemmi::Residue* res = nullptr;
        for (const gemmi::Residue& r : chain_a->residues) {
            if (r.seqid.num.value == 10) {
                res = &r;
                break;
            }
        }
        ASSERT_NE(res, nullptr);
        
        EXPECT_EQ(chain_a->name, "A");
        EXPECT_EQ(res->seqid.num.value, 10);
        EXPECT_EQ(res->name, "ILE");
        
        std::map<std::string, clipper::RTop_orth> rtops = coot::util::get_ori_to_this_res_gemmi(*res);
        EXPECT_EQ(rtops.size(), 1);
        
        std::map<std::string, clipper::RTop_orth>::const_iterator it = rtops.begin();
        EXPECT_EQ(it->first, "");  // empty alt conf
        // The oracle output showed rtop_valid: true, so we just verify the map is not empty
        EXPECT_FALSE(it->second.is_null());
    }

    // Case 2: Valid residue with mainchain atoms (GLY at A/11)
    {
        const gemmi::Model& model = st.models[0];
        const gemmi::Chain* chain_a = nullptr;
        for (const gemmi::Chain& chain : model.chains) {
            if (chain.name == "A") {
                chain_a = &chain;
                break;
            }
        }
        ASSERT_NE(chain_a, nullptr);
        
        const gemmi::Residue* res = nullptr;
        for (const gemmi::Residue& r : chain_a->residues) {
            if (r.seqid.num.value == 11) {
                res = &r;
                break;
            }
        }
        ASSERT_NE(res, nullptr);
        
        EXPECT_EQ(chain_a->name, "A");
        EXPECT_EQ(res->seqid.num.value, 11);
        EXPECT_EQ(res->name, "GLY");
        
        std::map<std::string, clipper::RTop_orth> rtops = coot::util::get_ori_to_this_res_gemmi(*res);
        EXPECT_EQ(rtops.size(), 1);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}