#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

using namespace coot;

TEST(OracleTest, crankshaft_set) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Case 1: Valid crankshaft_set with 4 consecutive residues
    // Using residues 10, 11, 12, 13 on chain A (ILE, GLY, GLU, GLY - all have N, C, O, CA)
    {
        gemmi::Model& model = st.models[0];
        gemmi::Chain& chain = model.chains[0];
        
        gemmi::Residue* res_0_ptr = nullptr;
        gemmi::Residue* res_1_ptr = nullptr;
        gemmi::Residue* res_2_ptr = nullptr;
        gemmi::Residue* res_3_ptr = nullptr;
        
        for (auto& res : chain.residues) {
            if (res.seqid.num.value == 10) res_0_ptr = &res;
            if (res.seqid.num.value == 11) res_1_ptr = &res;
            if (res.seqid.num.value == 12) res_2_ptr = &res;
            if (res.seqid.num.value == 13) res_3_ptr = &res;
        }
        
        ASSERT_NE(res_0_ptr, nullptr);
        ASSERT_NE(res_1_ptr, nullptr);
        ASSERT_NE(res_2_ptr, nullptr);
        ASSERT_NE(res_3_ptr, nullptr);

        EXPECT_NO_THROW({
            crankshaft_set cs(*res_0_ptr, *res_1_ptr, *res_2_ptr, *res_3_ptr);
            EXPECT_TRUE(cs.ca_1);
            EXPECT_TRUE(cs.ca_2);
            EXPECT_EQ(cs.v.size(), 8U);
        });
    }

    // Case 2: Invalid - try to use the same residue 4 times (should fail)
    {
        gemmi::Model& model = st.models[0];
        gemmi::Chain& chain = model.chains[0];
        gemmi::Residue& res_0 = chain.residues[9];
        
        EXPECT_THROW({
            crankshaft_set cs(res_0, res_0, res_0, res_0);
        }, std::runtime_error);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}