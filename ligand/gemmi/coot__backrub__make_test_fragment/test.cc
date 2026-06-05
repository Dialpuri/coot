#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, make_test_fragment) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model& model = st.models[0];

    // Find residues A/10, A/11, A/12
    const gemmi::Residue* prev_res = nullptr;
    const gemmi::Residue* this_res = nullptr;
    const gemmi::Residue* next_res = nullptr;

    for (gemmi::Chain& chain : model.chains) {
        if (chain.name != "A") continue;
        for (gemmi::Residue& res : chain.residues) {
            if (res.seqid.num.value == 10) prev_res = &res;
            if (res.seqid.num.value == 11) this_res = &res;
            if (res.seqid.num.value == 12) next_res = &res;
        }
    }

    ASSERT_NE(nullptr, prev_res);
    ASSERT_NE(nullptr, this_res);
    ASSERT_NE(nullptr, next_res);

    std::string chain_id_str = "A";

    // Find CA atoms in prev, this, next residues
    const gemmi::Atom* ca_prev_atom = nullptr;
    for (const gemmi::Atom& atm : prev_res->atoms) {
        std::string aname = atm.name;
        // Trim spaces for comparison
        if (aname == " CA" || aname == "CA" || aname == "ca") {
            ca_prev_atom = &atm;
            break;
        }
    }

    const gemmi::Atom* ca_next_atom = nullptr;
    for (const gemmi::Atom& atm : next_res->atoms) {
        std::string aname = atm.name;
        if (aname == " CA" || aname == "CA" || aname == "ca") {
            ca_next_atom = &atm;
            break;
        }
    }

    const gemmi::Atom* ca_this_atom = nullptr;
    for (const gemmi::Atom& atm : this_res->atoms) {
        std::string aname = atm.name;
        if (aname == " CA" || aname == "CA" || aname == "ca") {
            ca_this_atom = &atm;
            break;
        }
    }

    ASSERT_NE(nullptr, ca_prev_atom);
    ASSERT_NE(nullptr, ca_next_atom);
    ASSERT_NE(nullptr, ca_this_atom);

    clipper::Coord_orth ca_prev(ca_prev_atom->pos.x, ca_prev_atom->pos.y, ca_prev_atom->pos.z);
    clipper::Coord_orth ca_next(ca_next_atom->pos.x, ca_next_atom->pos.y, ca_next_atom->pos.z);

    // Case 1: 15 degree rotation
    {
        double rotation_angle = 15.0;
        coot::minimol::fragment frag = coot::backrub::make_test_fragment_gemmi(
            this_res, prev_res, next_res, chain_id_str,
            ca_prev, ca_next, rotation_angle);
        EXPECT_EQ(3, frag.n_filled_residues());
        EXPECT_EQ("A", frag.fragment_id);
        EXPECT_EQ(1, frag.min_res_no());
        EXPECT_EQ(12, frag.max_residue_number());
    }

    // Case 2: zero rotation (edge case)
    {
        double rotation_angle = 0.0;
        coot::minimol::fragment frag = coot::backrub::make_test_fragment_gemmi(
            this_res, prev_res, next_res, chain_id_str,
            ca_prev, ca_next, rotation_angle);
        EXPECT_EQ(3, frag.n_filled_residues());
        EXPECT_EQ("A", frag.fragment_id);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}