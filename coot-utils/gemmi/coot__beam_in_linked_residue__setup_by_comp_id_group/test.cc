#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, setup_by_comp_id_group) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Residue* asn_res = nullptr;

    for (auto& model : st.models) {
        for (auto& chain : model.chains) {
            for (auto& res : chain.residues) {
                if (res.name == "ASN") {
                    asn_res = &res;
                    break;
                }
            }
            if (asn_res) break;
        }
        if (asn_res) break;
    }

    if (!asn_res) {
        if (!st.models.empty()) {
            auto& model = st.models[0];
            if (!model.chains.empty()) {
                auto& chain = model.chains[0];
                if (!chain.residues.empty()) {
                    asn_res = &chain.residues[0];
                }
            }
        }
    }

    ASSERT_TRUE(asn_res != nullptr) << "No residues at all";

    std::string link_type = "glycosidic";
    bool have_template = false;
    gemmi::CRA template_res_ref{};
    gemmi::CRA template_res_mov{};

    // Case 1: comp_id_ref=ASN, group_new=N-sugar (file likely doesn't exist)
    {
        EXPECT_EQ(have_template, false);

        bool status = coot::beam_in_linked_residue::setup_by_comp_id_group_gemmi(
            "ASN", "N-sugar", link_type, template_res_ref, template_res_mov, have_template);

        EXPECT_EQ(status, false);
        EXPECT_EQ(have_template, false);
        EXPECT_EQ(template_res_ref.atom, nullptr);
        EXPECT_EQ(template_res_mov.atom, nullptr);
    }

    // Case 2: Invalid combination — file won't exist
    {
        link_type = "fake-link";
        have_template = false;
        template_res_ref = gemmi::CRA{};
        template_res_mov = gemmi::CRA{};

        bool status2 = coot::beam_in_linked_residue::setup_by_comp_id_group_gemmi(
            "ASN", "fake-group", link_type, template_res_ref, template_res_mov, have_template);

        EXPECT_EQ(status2, false);
        EXPECT_EQ(have_template, false);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}