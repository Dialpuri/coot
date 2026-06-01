#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, contact_atoms_info_t_contact_atom_t) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model& model = st.models[0];

    // Get two atoms: //A/1/CA (central) and //A/1/N (contactor) — same residue, different atoms
    gemmi::Residue* res1 = nullptr;
    for (auto& chain : model.chains) {
        if (chain.name == "A") {
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == 1) {
                    res1 = &res;
                    break;
                }
            }
            break;
        }
    }
    ASSERT_NE(res1, nullptr);

    gemmi::Atom* ca_atom = nullptr;
    gemmi::Atom* n_atom  = nullptr;
    for (auto& atom : res1->atoms) {
        if (atom.name == "CA") ca_atom = &atom;
        if (atom.name == "N") n_atom = &atom;
    }

    ASSERT_NE(ca_atom, nullptr);
    ASSERT_NE(n_atom, nullptr);

    // Build identity matrix
    double identity[4][4];
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            identity[i][j] = (i == j) ? 1.0 : 0.0;

    // Case 1: valid atoms — contactor = N, central = CA
    {
        auto con = coot::util::contact_atoms_info_t::contact_atom_t::contact_atom_t_gemmi(n_atom, ca_atom, identity);
        EXPECT_NEAR(con.dist, 1.4886, 1e-4);
        EXPECT_EQ(con.at->name, "N");
        EXPECT_NEAR(con.mat[0][0], 1.0, 1e-4);
        EXPECT_NEAR(con.mat[1][1], 1.0, 1e-4);
    }

    // Case 2: atoms further apart — use //A/1/CA and //A/10/CA
    {
        gemmi::Residue* res10 = nullptr;
        for (auto& chain : model.chains) {
            if (chain.name == "A") {
                for (auto& res : chain.residues) {
                    if (res.seqid.num.value == 10) {
                        res10 = &res;
                        break;
                    }
                }
                break;
            }
        }
        ASSERT_NE(res10, nullptr);

        gemmi::Atom* ca10 = nullptr;
        for (auto& atom : res10->atoms) {
            if (atom.name == "CA") ca10 = &atom;
        }
        ASSERT_NE(ca10, nullptr);

        auto con = coot::util::contact_atoms_info_t::contact_atom_t::contact_atom_t_gemmi(ca10, ca_atom, identity);
        EXPECT_NEAR(con.dist, 17.6110, 1e-4);
        EXPECT_EQ(con.at->name, "CA");
    }

    // Case 3: same atom as both — distance should be 0
    {
        auto con = coot::util::contact_atoms_info_t::contact_atom_t::contact_atom_t_gemmi(ca_atom, ca_atom, identity);
        EXPECT_NEAR(con.dist, 0.0, 1e-4);
        EXPECT_EQ(con.at->name, "CA");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}