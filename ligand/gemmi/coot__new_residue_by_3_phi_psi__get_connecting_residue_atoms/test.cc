#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, new_residue_by_3_phi_psi) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model& model = st.models[0];
    gemmi::Chain* chain = nullptr;
    for (auto& ch : model.chains) {
        if (ch.name == "A") { chain = &ch; break; }
    }
    ASSERT_NE(chain, nullptr);

    // Get residue 50 (actually A/59 ASN) as the main residue - has N, CA, C
    gemmi::Residue* res_50 = &chain->residues[49]; // 0-indexed
    gemmi::Residue* res_49 = &chain->residues[48]; // upstream (A/58 LEU)
    gemmi::Residue* res_51 = &chain->residues[50]; // downstream (A/60 HIS)

    // Case 1: residue with N, CA, C + upstream and downstream neighbours
    {
        coot::connecting_atoms_t_gemmi atoms =
            coot::get_connecting_residue_atoms_gemmi(res_50, res_49, res_51);

        EXPECT_FALSE(atoms.empty());
        EXPECT_NEAR(atoms.N_pos[0], 22.723, 1e-3);
        EXPECT_NEAR(atoms.N_pos[1], -8.957, 1e-3);
        EXPECT_NEAR(atoms.N_pos[2], 57.722, 1e-3);
        EXPECT_NEAR(atoms.CA_pos[0], 23.851, 1e-3);
        EXPECT_NEAR(atoms.CA_pos[1], -9.012, 1e-3);
        EXPECT_NEAR(atoms.CA_pos[2], 58.639, 1e-3);
        EXPECT_NEAR(atoms.C_pos[0], 25.133, 1e-3);
        EXPECT_NEAR(atoms.C_pos[1], -9.406, 1e-3);
        EXPECT_NEAR(atoms.C_pos[2], 57.918, 1e-3);

        EXPECT_TRUE(atoms.upstream_C.first);
        EXPECT_NEAR(atoms.upstream_C.second[0], 22.226, 1e-3);
        EXPECT_NEAR(atoms.upstream_C.second[1], -7.823, 1e-3);
        EXPECT_NEAR(atoms.upstream_C.second[2], 57.211, 1e-3);

        EXPECT_TRUE(atoms.downstream_N.first);
        EXPECT_NEAR(atoms.downstream_N.second[0], 26.184, 1e-3);
        EXPECT_NEAR(atoms.downstream_N.second[1], -8.589, 1e-3);
        EXPECT_NEAR(atoms.downstream_N.second[2], 58.081, 1e-3);

        auto phi = atoms.get_phi();
        auto psi = atoms.get_psi();
        EXPECT_TRUE(phi.first);
        EXPECT_NEAR(phi.second, -1.76229, 1e-4);
        EXPECT_TRUE(psi.first);
        EXPECT_NEAR(psi.second, 2.21248, 1e-4);
    }

    // Case 2: residue without upstream/downstream neighbours
    {
        coot::connecting_atoms_t_gemmi atoms =
            coot::get_connecting_residue_atoms_gemmi(res_50, nullptr, nullptr);

        EXPECT_FALSE(atoms.empty());
        EXPECT_NEAR(atoms.N_pos[0], 22.723, 1e-3);
        EXPECT_NEAR(atoms.N_pos[1], -8.957, 1e-3);
        EXPECT_NEAR(atoms.N_pos[2], 57.722, 1e-3);
        EXPECT_FALSE(atoms.upstream_C.first);
        EXPECT_FALSE(atoms.downstream_N.first);
    }

    // Case 3: Glycine residue (A/11 GLY) - has N, CA, C but no CB
    {
        gemmi::Residue* res_gly = &chain->residues[10]; // GLY 11
        gemmi::Residue* res_gly_up = &chain->residues[9]; // ILE 10
        gemmi::Residue* res_gly_down = &chain->residues[11]; // GLU 12

        coot::connecting_atoms_t_gemmi atoms =
            coot::get_connecting_residue_atoms_gemmi(res_gly, res_gly_up, res_gly_down);

        EXPECT_FALSE(atoms.empty());
        EXPECT_NEAR(atoms.N_pos[0], 25.736, 1e-3);
        EXPECT_NEAR(atoms.N_pos[1], 13.244, 1e-3);
        EXPECT_NEAR(atoms.N_pos[2], 66.003, 1e-3);
        EXPECT_NEAR(atoms.CA_pos[0], 25.487, 1e-3);
        EXPECT_NEAR(atoms.CA_pos[1], 13.854, 1e-3);
        EXPECT_NEAR(atoms.CA_pos[2], 64.707, 1e-3);
        EXPECT_NEAR(atoms.C_pos[0], 24.004, 1e-3);
        EXPECT_NEAR(atoms.C_pos[1], 13.943, 1e-3);
        EXPECT_NEAR(atoms.C_pos[2], 64.377, 1e-3);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}