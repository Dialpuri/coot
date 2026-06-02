#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, dict_link_info_t) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    coot::protein_geometry geom;
    geom.init_standard();

    // Helper: find chain and residue by chain name and sequence number
    auto find_residue = [&](const std::string &chain_name, int seqnum) -> std::pair<const gemmi::Chain*, const gemmi::Residue*> {
        for (gemmi::Model &model : st.models) {
            for (gemmi::Chain &chain : model.chains) {
                if (chain.name == chain_name) {
                    for (gemmi::Residue &res : chain.residues) {
                        if (res.seqid.num.value == seqnum) {
                            return {&chain, &res};
                        }
                    }
                }
            }
        }
        return {nullptr, nullptr};
    };

    // Case 1: Valid residues but link not found in dictionary
    {
        auto cr = find_residue("A", 5);
        const gemmi::Chain *chain_ref = cr.first;
        const gemmi::Residue *res_ref = cr.second;
        auto cn = find_residue("A", 6);
        const gemmi::Chain *chain_new = cn.first;
        const gemmi::Residue *res_new = cn.second;

        EXPECT_NE(res_ref, nullptr);
        EXPECT_NE(res_new, nullptr);

        bool got_exception = false;
        std::string msg;
        try {
            coot::dict_link_info_t info = coot::dict_link_info_t::dict_link_info_t_gemmi(chain_ref, res_ref, chain_new, res_new, "PEPT", geom);
        } catch (const std::runtime_error &e) {
            got_exception = true;
            msg = e.what();
        }
        EXPECT_TRUE(got_exception);
        EXPECT_NE(msg.find("Link not found"), std::string::npos);
        EXPECT_NE(msg.find("PEPT"), std::string::npos);
        EXPECT_NE(msg.find("GLN"), std::string::npos);
        EXPECT_NE(msg.find("LYS"), std::string::npos);
    }

    // Case 2: Non-existent residue — res_new is null
    {
        auto cr = find_residue("A", 5);
        const gemmi::Residue *res_ref = cr.second;
        auto cn = find_residue("A", 9999);
        const gemmi::Residue *res_new = cn.second;

        EXPECT_NE(res_ref, nullptr);
        EXPECT_EQ(res_new, nullptr);
    }

    // Case 3: Null residue_ref — throws "Null residue_ref"
    {
        const gemmi::Chain *chain_ref = nullptr;
        const gemmi::Residue *res_ref = nullptr;
        auto cn = find_residue("A", 6);
        const gemmi::Chain *chain_new = cn.first;
        const gemmi::Residue *res_new = cn.second;

        EXPECT_EQ(res_ref, nullptr);
        EXPECT_NE(res_new, nullptr);

        bool got_exception = false;
        std::string msg;
        try {
            coot::dict_link_info_t info = coot::dict_link_info_t::dict_link_info_t_gemmi(chain_ref, res_ref, chain_new, res_new, "PEPT", geom);
        } catch (const std::runtime_error &e) {
            got_exception = true;
            msg = e.what();
        }
        EXPECT_TRUE(got_exception);
        EXPECT_NE(msg.find("Null residue_ref"), std::string::npos);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}