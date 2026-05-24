#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, IsSSBonded) {
    // Load the PDB file
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    const gemmi::Model& model = st.models[0];
    
    // Test case 1: CYS residue 118 in chain A
    {
        const gemmi::Residue* res118 = nullptr;
        std::string chain_id_118;
        for (const auto& chain : model.chains) {
            for (const auto& res : chain.residues) {
                if (res.seqid.num.value == 118 && chain.name == "A") {
                    res118 = &res;
                    chain_id_118 = chain.name;
                    break;
                }
            }
            if (res118) break;
        }
        ASSERT_NE(res118, nullptr);
        bool is_ss = coot::reduce::is_ss_bonded_gemmi(res118, model);
        EXPECT_FALSE(is_ss);
        EXPECT_EQ(chain_id_118, "A");
        EXPECT_EQ(res118->seqid.num.value, 118);
        EXPECT_EQ(res118->name, "CYS");
    }
    
    // Test case 2: CYS residue 177 in chain A
    {
        const gemmi::Residue* res177 = nullptr;
        std::string chain_id_177;
        for (const auto& chain : model.chains) {
            for (const auto& res : chain.residues) {
                if (res.seqid.num.value == 177 && chain.name == "A") {
                    res177 = &res;
                    chain_id_177 = chain.name;
                    break;
                }
            }
            if (res177) break;
        }
        ASSERT_NE(res177, nullptr);
        bool is_ss = coot::reduce::is_ss_bonded_gemmi(res177, model);
        EXPECT_FALSE(is_ss);
        EXPECT_EQ(chain_id_177, "A");
        EXPECT_EQ(res177->seqid.num.value, 177);
        EXPECT_EQ(res177->name, "CYS");
    }
    
    // Test case 3: Non-CYS residue 10 in chain A
    {
        const gemmi::Residue* res10 = nullptr;
        std::string chain_id_10;
        for (const auto& chain : model.chains) {
            for (const auto& res : chain.residues) {
                if (res.seqid.num.value == 10 && chain.name == "A") {
                    res10 = &res;
                    chain_id_10 = chain.name;
                    break;
                }
            }
            if (res10) break;
        }
        ASSERT_NE(res10, nullptr);
        bool is_ss = coot::reduce::is_ss_bonded_gemmi(res10, model);
        EXPECT_FALSE(is_ss);
        EXPECT_EQ(chain_id_10, "A");
        EXPECT_EQ(res10->seqid.num.value, 10);
        EXPECT_EQ(res10->name, "ILE");
    }
    
    // Test case 4: NULL residue (edge case)
    {
        const gemmi::Residue* null_res = nullptr;
        bool is_ss = coot::reduce::is_ss_bonded_gemmi(null_res, model);
        EXPECT_FALSE(is_ss);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}