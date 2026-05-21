#define HAVE_BOOST_BASED_THREAD_POOL_LIBRARY

#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, peptide_C_and_N_are_in_order_p) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];
    
    // case: Two consecutive residues in proper order (GLU 2 -> ASN 3)
    {
        gemmi::Residue* res1 = nullptr;
        gemmi::Residue* res2 = nullptr;
        for (auto& chain : model.chains) {
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == 2 && res.name == "GLU") {
                    res1 = &res;
                }
                if (res.seqid.num.value == 3 && res.name == "ASN") {
                    res2 = &res;
                }
            }
        }
        ASSERT_NE(res1, nullptr);
        ASSERT_NE(res2, nullptr);
        
        gemmi::CRA cra1{&model.chains[0], res1, nullptr};
        gemmi::CRA cra2{&model.chains[0], res2, nullptr};
        
        auto result = coot::peptide_C_and_N_are_in_order_p_gemmi(cra1, cra2);
        auto info = result.first;
        bool switch_needed = result.second;
        
        EXPECT_EQ(res1->name, "GLU");
        EXPECT_EQ(res2->name, "ASN");
        EXPECT_EQ(info, coot::restraints_container_t::IS_PEPTIDE);
        EXPECT_FALSE(switch_needed);
    }
    
    // case: Two consecutive residues where res2 comes before res1 in serial order (ASN 3 -> GLU 2)
    {
        gemmi::Residue* res1 = nullptr;
        gemmi::Residue* res2 = nullptr;
        for (auto& chain : model.chains) {
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == 3 && res.name == "ASN") {
                    res1 = &res;
                }
                if (res.seqid.num.value == 2 && res.name == "GLU") {
                    res2 = &res;
                }
            }
        }
        ASSERT_NE(res1, nullptr);
        ASSERT_NE(res2, nullptr);
        
        gemmi::CRA cra1{&model.chains[0], res1, nullptr};
        gemmi::CRA cra2{&model.chains[0], res2, nullptr};
        
        auto result = coot::peptide_C_and_N_are_in_order_p_gemmi(cra1, cra2);
        auto info = result.first;
        bool switch_needed = result.second;
        
        EXPECT_EQ(res1->name, "ASN");
        EXPECT_EQ(res2->name, "GLU");
        EXPECT_EQ(info, coot::restraints_container_t::IS_PEPTIDE);
        EXPECT_TRUE(switch_needed);
    }
    
    // case: Non-consecutive residues (GLU 2 -> GLN 5)
    {
        gemmi::Residue* res1 = nullptr;
        gemmi::Residue* res2 = nullptr;
        for (auto& chain : model.chains) {
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == 2 && res.name == "GLU") {
                    res1 = &res;
                }
                if (res.seqid.num.value == 5 && res.name == "GLN") {
                    res2 = &res;
                }
            }
        }
        ASSERT_NE(res1, nullptr);
        ASSERT_NE(res2, nullptr);
        
        gemmi::CRA cra1{&model.chains[0], res1, nullptr};
        gemmi::CRA cra2{&model.chains[0], res2, nullptr};
        
        auto result = coot::peptide_C_and_N_are_in_order_p_gemmi(cra1, cra2);
        auto info = result.first;
        bool switch_needed = result.second;
        
        EXPECT_EQ(res1->name, "GLU");
        EXPECT_EQ(res2->name, "GLN");
        EXPECT_EQ(info, coot::restraints_container_t::IS_NOT_PEPTIDE);
        EXPECT_FALSE(switch_needed);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}