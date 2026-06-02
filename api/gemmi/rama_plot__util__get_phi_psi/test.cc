#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_phi_psi) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model& model = st.models[0];
    ASSERT_GE(model.chains.size(), 1u);

    gemmi::Chain& chain = model.chains[0];

    size_t nres = chain.residues.size();
    EXPECT_EQ(nres, 267u);

    // Case 1: typical valid residue in the middle (indices 1,2,3)
    {
        gemmi::Residue* res_prev = &chain.residues[1];
        gemmi::Residue* res_this = &chain.residues[2];
        gemmi::Residue* res_next = &chain.residues[3];

        ASSERT_NE(res_prev, nullptr);
        ASSERT_NE(res_this, nullptr);
        ASSERT_NE(res_next, nullptr);

        EXPECT_EQ(res_prev->seqid.num.value, 2);
        EXPECT_STREQ(res_prev->name.c_str(), "GLU");
        EXPECT_EQ(res_this->seqid.num.value, 3);
        EXPECT_STREQ(res_this->name.c_str(), "ASN");
        EXPECT_EQ(res_next->seqid.num.value, 4);
        EXPECT_STREQ(res_next->name.c_str(), "PHE");

        std::pair<bool, rama_plot::phi_psi_t> result =
            rama_plot::util::get_phi_psi_gemmi(
                gemmi::CRA{&chain, res_prev, nullptr},
                gemmi::CRA{&chain, res_this, nullptr},
                gemmi::CRA{&chain, res_next, nullptr});

        EXPECT_TRUE(result.first);
        EXPECT_NEAR(result.second.phi, -75.3333, 1e-4);
        EXPECT_NEAR(result.second.psi, 1.75129, 1e-4);
        EXPECT_EQ(result.second.label, "3 A ASN");
        EXPECT_FALSE(result.second.is_pre_pro);
    }

    // Case 2: boundary case (indices 0,1,2)
    {
        gemmi::Residue* res_prev = &chain.residues[0];
        gemmi::Residue* res_this = &chain.residues[1];
        gemmi::Residue* res_next = &chain.residues[2];

        ASSERT_NE(res_prev, nullptr);
        ASSERT_NE(res_this, nullptr);
        ASSERT_NE(res_next, nullptr);

        EXPECT_EQ(res_prev->seqid.num.value, 1);
        EXPECT_STREQ(res_prev->name.c_str(), "MET");
        EXPECT_EQ(res_this->seqid.num.value, 2);
        EXPECT_STREQ(res_this->name.c_str(), "GLU");
        EXPECT_EQ(res_next->seqid.num.value, 3);
        EXPECT_STREQ(res_next->name.c_str(), "ASN");

        std::pair<bool, rama_plot::phi_psi_t> result =
            rama_plot::util::get_phi_psi_gemmi(
                gemmi::CRA{&chain, res_prev, nullptr},
                gemmi::CRA{&chain, res_this, nullptr},
                gemmi::CRA{&chain, res_next, nullptr});

        EXPECT_TRUE(result.first);
        EXPECT_NEAR(result.second.phi, -55.9405, 1e-4);
        EXPECT_NEAR(result.second.psi, -13.5155, 1e-4);
    }

    // Case 3: PRO residue (proline has special handling)
    {
        gemmi::Residue* pro_prev = nullptr;
        gemmi::Residue* pro_this = nullptr;
        gemmi::Residue* pro_next = nullptr;

        for (size_t i = 1; i + 1 < nres; i++) {
            gemmi::Residue* r = &chain.residues[i];
            if (r && std::string(r->name) == "PRO") {
                pro_prev = &chain.residues[i - 1];
                pro_this = r;
                pro_next = &chain.residues[i + 1];
                if (pro_prev && pro_next) break;
            }
        }

        ASSERT_NE(pro_this, nullptr);

        EXPECT_EQ(pro_prev->seqid.num.value, 60);
        EXPECT_STREQ(pro_prev->name.c_str(), "HIS");
        EXPECT_EQ(pro_this->seqid.num.value, 61);
        EXPECT_STREQ(pro_this->name.c_str(), "PRO");
        EXPECT_EQ(pro_next->seqid.num.value, 62);
        EXPECT_STREQ(pro_next->name.c_str(), "ASN");

        std::pair<bool, rama_plot::phi_psi_t> result =
            rama_plot::util::get_phi_psi_gemmi(
                gemmi::CRA{&chain, pro_prev, nullptr},
                gemmi::CRA{&chain, pro_this, nullptr},
                gemmi::CRA{&chain, pro_next, nullptr});

        EXPECT_TRUE(result.first);
        EXPECT_NEAR(result.second.phi, -58.962, 1e-4);
        EXPECT_NEAR(result.second.psi, -20.2551, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}