#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, transform_chain) {
    // Setup: load PDB and get chain A + first CA atom
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty()) << "no models";

    gemmi::Model& model = st.models[0];
    ASSERT_FALSE(model.chains.empty()) << "no chains";

    gemmi::Chain& chain_A = model.chains[0];
    ASSERT_FALSE(chain_A.residues.empty()) << "no residues in chain A";

    // Get first CA atom to track position changes
    gemmi::Residue& first_res = chain_A.residues[0];
    const gemmi::Atom* ca_atom = nullptr;
    for (const gemmi::Atom& a : first_res.atoms) {
        if (a.name.find("CA") != std::string::npos) {
            ca_atom = &a;
            break;
        }
    }
    ASSERT_NE(ca_atom, nullptr) << "null CA atom in first residue";

    // Case 1: Identity transformation (no change expected)
    {
        double x_before = ca_atom->pos.x;
        double y_before = ca_atom->pos.y;
        double z_before = ca_atom->pos.z;

        clipper::Mat33<double> identity(1, 0, 0, 0, 1, 0, 0, 0, 1);
        clipper::Coord_orth zero_trans(0.0, 0.0, 0.0);
        clipper::RTop_orth identity_rtop(identity, zero_trans);

        coot::util::transform_chain_gemmi(chain_A, identity_rtop);

        EXPECT_NEAR(ca_atom->pos.x, 13.401, 1e-4);
        EXPECT_NEAR(ca_atom->pos.y, 3.78, 1e-4);
        EXPECT_NEAR(ca_atom->pos.z, 74.703, 1e-4);
    }

    // Case 2: Translation by (5, 0, 0) — positions should change
    {
        double x2_before = ca_atom->pos.x;
        double y2_before = ca_atom->pos.y;
        double z2_before = ca_atom->pos.z;

        EXPECT_NEAR(x2_before, 13.401, 1e-4);
        EXPECT_NEAR(y2_before, 3.78, 1e-4);
        EXPECT_NEAR(z2_before, 74.703, 1e-4);

        clipper::Mat33<double> identity(1, 0, 0, 0, 1, 0, 0, 0, 1);
        clipper::Coord_orth trans(5.0, 0.0, 0.0);
        clipper::RTop_orth trans_rtop(identity, trans);

        coot::util::transform_chain_gemmi(chain_A, trans_rtop);

        EXPECT_NEAR(ca_atom->pos.x, 18.401, 1e-4);
        EXPECT_NEAR(ca_atom->pos.y, 3.78, 1e-4);
        EXPECT_NEAR(ca_atom->pos.z, 74.703, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}