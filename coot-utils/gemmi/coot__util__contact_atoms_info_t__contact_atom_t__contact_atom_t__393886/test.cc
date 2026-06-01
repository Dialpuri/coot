#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

static gemmi::Atom* find_atom(const gemmi::Model& model,
                              const std::string& chain_id,
                              int res_seq_num,
                              const std::string& atom_name) {
  for (const gemmi::Chain& ch : model.chains) {
    if (ch.name != chain_id) continue;
    for (const gemmi::Residue& res : ch.residues) {
      if (res.seqid.num.value != res_seq_num) continue;
      for (const gemmi::Atom& atm : res.atoms) {
        if (atm.name == atom_name) return const_cast<gemmi::Atom*>(&atm);
      }
    }
  }
  return nullptr;
}

TEST(OracleTest, contact_atoms_info_t_contact_atom_t) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_EQ(st.models.size(), 1u);
    const gemmi::Model &model = st.models[0];

    // Case 1: Two different atoms — normal case
    {
        gemmi::Atom *a1 = find_atom(model, "A", 1, "CA");
        gemmi::Atom *a2 = find_atom(model, "A", 1, "N");
        ASSERT_NE(a1, nullptr);
        ASSERT_NE(a2, nullptr);

        EXPECT_NEAR(a1->pos.x, 13.401, 1e-4);
        EXPECT_NEAR(a1->pos.y, 3.78, 1e-4);
        EXPECT_NEAR(a1->pos.z, 74.703, 1e-4);
        EXPECT_NEAR(a2->pos.x, 12.334, 1e-4);
        EXPECT_NEAR(a2->pos.y, 2.772, 1e-4);
        EXPECT_NEAR(a2->pos.z, 74.951, 1e-4);

        coot::util::contact_atoms_info_t::contact_atom_t ca(a1, a2);
        EXPECT_NEAR(ca.dist, 1.48864, 1e-4);
        // Identity matrix on diagonal
        EXPECT_NEAR(ca.mat[0][0], 1, 1e-4);
        EXPECT_NEAR(ca.mat[1][1], 1, 1e-4);
        EXPECT_NEAR(ca.mat[2][2], 1, 1e-4);
        EXPECT_NEAR(ca.mat[3][3], 1, 1e-4);
    }

    // Case 2: Same atom for both — distance should be 0
    {
        gemmi::Atom *a = find_atom(model, "A", 1, "CA");
        ASSERT_NE(a, nullptr);

        coot::util::contact_atoms_info_t::contact_atom_t ca(a, a);
        EXPECT_NEAR(ca.dist, 0, 1e-4);
    }

    // Case 3: With matrix overload — pass a non-identity matrix
    {
        gemmi::Atom *a1 = find_atom(model, "A", 2, "CA");
        gemmi::Atom *a2 = find_atom(model, "A", 3, "CA");
        ASSERT_NE(a1, nullptr);
        ASSERT_NE(a2, nullptr);

        double mat_in[4][4];
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                mat_in[i][j] = (i == j) ? 2.0 : 0.5;

        EXPECT_NEAR(mat_in[0][0], 2, 1e-4);
        EXPECT_NEAR(mat_in[1][1], 2, 1e-4);

        coot::util::contact_atoms_info_t::contact_atom_t ca(a1, a2, mat_in);
        EXPECT_NEAR(ca.dist, 3.81995, 1e-4);
        EXPECT_NEAR(ca.mat[0][0], 2, 1e-4);
        EXPECT_NEAR(ca.mat[1][1], 2, 1e-4);
        EXPECT_NEAR(ca.mat[0][1], 0.5, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}