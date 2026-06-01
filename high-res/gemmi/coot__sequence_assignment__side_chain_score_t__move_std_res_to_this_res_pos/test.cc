#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include "clipper/core/clipper_types.h"
#include "clipper/core/coords.h"

TEST(OracleTest, move_std_res_to_this_res_pos) {
    // Load PDB and get residue from chain A
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model *model = &st.models[0];
    ASSERT_NE(model, nullptr);

    gemmi::Chain *chain = nullptr;
    for (gemmi::Chain &ch : model->chains) {
        if (ch.name == "A") {
            chain = &ch;
            break;
        }
    }
    ASSERT_NE(chain, nullptr);

    gemmi::Residue *res = &chain->residues[0];
    ASSERT_NE(res, nullptr);

    EXPECT_EQ(std::string(res->name), "MET");
    EXPECT_EQ(res->seqid.num.value, 1);

    int nAtoms = static_cast<int>(res->atoms.size());
    EXPECT_EQ(nAtoms, 8);

    // Record coordinates of first atom BEFORE identity transform
    double x_before = res->atoms[0].pos.x;
    double y_before = res->atoms[0].pos.y;
    double z_before = res->atoms[0].pos.z;

    // gemmi atom.name is unpadded ("N"); MMDB returned " N  " (4-char padded)
    // PDB atom names are right-justified in 4 chars; single-char element names
    // get a leading space from gemmi::padded_name().
    std::string mmdb_style_name = res->atoms[0].padded_name();
    while (mmdb_style_name.size() < 4) mmdb_style_name += ' ';
    EXPECT_EQ(mmdb_style_name, " N  ");
    EXPECT_NEAR(x_before, 12.334, 1e-4);
    EXPECT_NEAR(y_before, 2.772, 1e-4);
    EXPECT_NEAR(z_before, 74.951, 1e-4);

    // Create side_chain_score_t instance
    coot::sequence_assignment::side_chain_score_t scorer;

    // Case 1: Identity transform (no change expected)
    {
        clipper::Mat33<double> id_mat(1, 0, 0, 0, 1, 0, 0, 0, 1);
        clipper::Coord_orth zero_trn(0, 0, 0);
        clipper::RTop_orth identity_rtop(id_mat, zero_trn);
        EXPECT_NO_THROW(scorer.move_std_res_to_this_res_pos_gemmi(identity_rtop, res));

        EXPECT_NEAR(res->atoms[0].pos.x, 12.334, 1e-4);
        EXPECT_NEAR(res->atoms[0].pos.y, 2.772, 1e-4);
        EXPECT_NEAR(res->atoms[0].pos.z, 74.951, 1e-4);
    }

    // Case 2: Translation by (10, 20, 30) — coordinates should change
    {
        clipper::Mat33<double> id_mat(1, 0, 0, 0, 1, 0, 0, 0, 1);
        clipper::Coord_orth trans_vec(10.0, 20.0, 30.0);
        clipper::RTop_orth trans_rtop(id_mat, trans_vec);
        EXPECT_NO_THROW(scorer.move_std_res_to_this_res_pos_gemmi(trans_rtop, res));

        EXPECT_NEAR(res->atoms[0].pos.x, 22.334, 1e-4);
        EXPECT_NEAR(res->atoms[0].pos.y, 22.772, 1e-4);
        EXPECT_NEAR(res->atoms[0].pos.z, 104.951, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}