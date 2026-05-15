#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, co_gemmi) {
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // case 1: valid atom — chain A, residue seqnum 10, atom name CA
    {
        gemmi::Atom *atom = nullptr;
        for (auto &model : st.models) {
            for (auto &chain : model.chains) {
                if (chain.name != "A") continue;
                for (auto &res : chain.residues) {
                    if (res.seqid.num.value != 10) continue;
                    for (auto &atm : res.atoms) {
                        if (atm.name == "CA") {
                            atom = &atm;
                            goto found;
                        }
                    }
                }
            }
        }
    found:
        ASSERT_NE(atom, nullptr);
        clipper::Coord_orth pt = coot::lidia_utils::co_gemmi(atom);
        EXPECT_NEAR(pt.x(), 27.175, 1e-4);
        EXPECT_NEAR(pt.y(), 12.237, 1e-4);
        EXPECT_NEAR(pt.z(), 67.71, 1e-4);
    }

    // case 2: invalid CID — residue 9999 does not exist
    {
        gemmi::Atom *atom = nullptr;
        for (auto &model : st.models) {
            for (auto &chain : model.chains) {
                if (chain.name != "A") continue;
                for (auto &res : chain.residues) {
                    if (res.seqid.num.value != 9999) continue;
                    for (auto &atm : res.atoms) {
                        if (atm.name == "N") {
                            atom = &atm;
                            goto found2;
                        }
                    }
                }
            }
        }
    found2:
        EXPECT_EQ(atom, nullptr);
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}