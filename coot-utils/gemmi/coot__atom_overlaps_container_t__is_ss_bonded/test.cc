#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, is_ss_bonded) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Find a CYS residue and a non-CYS residue (track chain name for assertions)
    const gemmi::Residue* cys_res = nullptr;
    const gemmi::Residue* non_cys_res = nullptr;
    std::string cys_chain_name, non_cys_chain_name;

    if (st.models.empty()) {
        FAIL() << "No models";
    }

    gemmi::Model& model = st.models[0];

    for (gemmi::Chain& chain : model.chains) {
        for (gemmi::Residue& res : chain.residues) {
            std::string rname = res.name;
            if (rname == "CYS" && cys_res == nullptr) {
                cys_res = &res;
                cys_chain_name = chain.name;
            }
            if (rname != "CYS" && non_cys_res == nullptr) {
                non_cys_res = &res;
                non_cys_chain_name = chain.name;
            }
        }
    }

    if (non_cys_res == nullptr) {
        FAIL() << "No non-CYS residue found";
    }

    // Case 1: CYS residue (should return true)
    {
        if (cys_res) {
            std::string res_name = cys_res->name;
            int seqnum = cys_res->seqid.num.value;
            EXPECT_EQ(cys_chain_name, "A");
            EXPECT_EQ(seqnum, 118);
            EXPECT_EQ(res_name, "CYS");
            bool result = coot::is_ss_bonded_gemmi(cys_res);
            EXPECT_TRUE(result);
        }
    }

    // Case 2: non-CYS residue (should return false)
    {
        std::string res_name = non_cys_res->name;
        int seqnum = non_cys_res->seqid.num.value;
        EXPECT_EQ(non_cys_chain_name, "A");
        EXPECT_EQ(seqnum, 1);
        EXPECT_EQ(res_name, "MET");
        bool result = coot::is_ss_bonded_gemmi(non_cys_res);
        EXPECT_FALSE(result);
    }

    // Case 3: nullptr residue (should return false — guarded path)
    {
        bool result = coot::is_ss_bonded_gemmi(nullptr);
        EXPECT_FALSE(result);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}