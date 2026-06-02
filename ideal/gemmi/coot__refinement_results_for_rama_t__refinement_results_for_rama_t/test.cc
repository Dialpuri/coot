#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

static const gemmi::Atom *find_atom_by_cid(gemmi::Model &model, const std::string &cid) {
    // CID format: "//<chain>/<resno>/<atomname>"
    if (cid.size() < 3 || cid.substr(0, 2) != "//") return nullptr;
    std::string rest = cid.substr(2);
    size_t pos1 = rest.find('/');
    if (pos1 == std::string::npos) return nullptr;
    size_t pos2 = rest.find('/', pos1 + 1);
    if (pos2 == std::string::npos) return nullptr;
    std::string chain_id = rest.substr(0, pos1);
    std::string resno_str = rest.substr(pos1 + 1, pos2 - pos1 - 1);
    std::string atom_name = rest.substr(pos2 + 1);
    int resno = std::stoi(resno_str);

    for (auto &chain : model.chains) {
        if (chain.name != chain_id) continue;
        for (auto &res : chain.residues) {
            if (res.seqid.num.value != resno) continue;
            for (auto &atom : res.atoms) {
                if (atom.name == atom_name) {
                    return &atom;
                }
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, refinement_results_for_rama) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model &model = st.models[0];

    // Case 1: All 5 atoms present (typical backbone: C(i-1), N(i), CA(i), C(i), N(i+1))
    {
        std::string cid1 = "//A/10/C";
        std::string cid2 = "//A/11/N";
        std::string cid3 = "//A/11/CA";
        std::string cid4 = "//A/11/C";
        std::string cid5 = "//A/12/N";

        const gemmi::Atom *at1 = find_atom_by_cid(model, cid1);
        const gemmi::Atom *at2 = find_atom_by_cid(model, cid2);
        const gemmi::Atom *at3 = find_atom_by_cid(model, cid3);
        const gemmi::Atom *at4 = find_atom_by_cid(model, cid4);
        const gemmi::Atom *at5 = find_atom_by_cid(model, cid5);

        EXPECT_NE(at1, nullptr);
        EXPECT_NE(at2, nullptr);
        EXPECT_NE(at3, nullptr);
        EXPECT_NE(at4, nullptr);
        EXPECT_NE(at5, nullptr);

        float distortion_in = 12.5f;

        coot::refinement_results_for_rama_t_gemmi rr(model, at1, at2, at3, at4, at5, distortion_in);

        EXPECT_NEAR(rr.distortion, 12.5f, 1e-4);
        EXPECT_EQ(rr.atom_spec_CA.format(), "[spec: model 1 \"A\"   11 \"\" \" CA \" \"\"]");
        EXPECT_NEAR(rr.ball_pos_x, 25.5949f, 1e-3);
        EXPECT_NEAR(rr.ball_pos_y, 14.2387f, 1e-3);
        EXPECT_NEAR(rr.ball_pos_z, 64.7699f, 1e-3);
    }

    // Case 2: Only at_3 (CA) present — tests partial guard path
    {
        std::string cid3 = "//A/11/CA";
        const gemmi::Atom *at3 = find_atom_by_cid(model, cid3);

        EXPECT_NE(at3, nullptr);

        float distortion_in = 5.0f;

        coot::refinement_results_for_rama_t_gemmi rr(model, nullptr, nullptr, at3, nullptr, nullptr, distortion_in);

        EXPECT_NEAR(rr.distortion, 5.0f, 1e-4);
        EXPECT_EQ(rr.atom_spec_CA.format(), "[spec: model 1 \"A\"   11 \"\" \" CA \" \"\"]");
        EXPECT_NEAR(rr.ball_pos_x, 25.987f, 1e-3);
        EXPECT_NEAR(rr.ball_pos_y, 13.854f, 1e-3);
        EXPECT_NEAR(rr.ball_pos_z, 64.707f, 1e-3);
    }

    // Case 3: All atoms null — ball_pos stays at (0,0,0)
    {
        float distortion_in = 0.0f;

        coot::refinement_results_for_rama_t_gemmi rr(model, nullptr, nullptr, nullptr, nullptr, nullptr, distortion_in);

        EXPECT_NEAR(rr.distortion, 0.0f, 1e-4);
        EXPECT_EQ(rr.atom_spec_CA.format(), "[spec: model -1 \"unset\" -2147483647 \"\" \"\" \"\"]");
        EXPECT_NEAR(rr.ball_pos_x, 0.0f, 1e-4);
        EXPECT_NEAR(rr.ball_pos_y, 0.0f, 1e-4);
        EXPECT_NEAR(rr.ball_pos_z, 0.0f, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}