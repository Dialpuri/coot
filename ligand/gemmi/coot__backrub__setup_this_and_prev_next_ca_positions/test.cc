#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, backrub) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model& model = st.models[0];

    // Find chain A
    const gemmi::Chain* chain_a = nullptr;
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == "A") {
            chain_a = &chain;
            break;
        }
    }
    ASSERT_NE(chain_a, nullptr);

    // Find residues 2, 3, 4
    const gemmi::Residue* prev_res = nullptr;
    const gemmi::Residue* this_res = nullptr;
    const gemmi::Residue* next_res = nullptr;

    for (const gemmi::Residue& res : chain_a->residues) {
        if (res.seqid.num.value == 2) prev_res = &res;
        if (res.seqid.num.value == 3) this_res = &res;
        if (res.seqid.num.value == 4) next_res = &res;
    }

    EXPECT_NE(prev_res, nullptr);
    EXPECT_NE(this_res, nullptr);
    EXPECT_NE(next_res, nullptr);

    // Case 1: Valid three consecutive residues (prev=2/GLU, this=3/ASN, next=4/PHE)
    {
        clipper::Coord_orth ca_this(0, 0, 0);
        clipper::Coord_orth ca_prev(0, 0, 0);
        clipper::Coord_orth ca_next(0, 0, 0);

        coot::backrub::setup_this_and_prev_next_ca_positions_gemmi(
            this_res, prev_res, next_res, "", ca_this, ca_prev, ca_next);

        EXPECT_NEAR(ca_this.x(), 15.902, 1e-4);
        EXPECT_NEAR(ca_this.y(), 1.647, 1e-4);
        EXPECT_NEAR(ca_this.z(), 78.426, 1e-4);
        EXPECT_NEAR(ca_prev.x(), 13.55, 1e-4);
        EXPECT_NEAR(ca_prev.y(), 4.657, 1e-4);
        EXPECT_NEAR(ca_prev.z(), 78.422, 1e-4);
        EXPECT_NEAR(ca_next.x(), 18.723, 1e-4);
        EXPECT_NEAR(ca_next.y(), 3.91, 1e-4);
        EXPECT_NEAR(ca_next.z(), 77.163, 1e-4);
    }

    // Case 2: Null prev residue — call with nullptr
    {
        clipper::Coord_orth ca_this(0, 0, 0);
        clipper::Coord_orth ca_prev(0, 0, 0);
        clipper::Coord_orth ca_next(0, 0, 0);

        EXPECT_THROW(
            coot::backrub::setup_this_and_prev_next_ca_positions_gemmi(
                this_res, nullptr, next_res, "", ca_this, ca_prev, ca_next),
            std::runtime_error);
    }

    // Case 3: Null next residue
    {
        clipper::Coord_orth ca_this(0, 0, 0);
        clipper::Coord_orth ca_prev(0, 0, 0);
        clipper::Coord_orth ca_next(0, 0, 0);

        EXPECT_THROW(
            coot::backrub::setup_this_and_prev_next_ca_positions_gemmi(
                this_res, prev_res, nullptr, "", ca_this, ca_prev, ca_next),
            std::runtime_error);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}