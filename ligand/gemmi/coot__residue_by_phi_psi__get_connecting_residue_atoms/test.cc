#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__molecule_t__get_residue/gemmi/function.hh"

TEST(OracleTest, residue_by_phi_psi) {

    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];

    // Case 1: Middle residue with both neighbours
    {
        gemmi::Residue *res = coot::get_residue_gemmi("//A/100", model);
        gemmi::Residue *upstream = coot::get_residue_gemmi("//A/99", model);
        gemmi::Residue *downstream = coot::get_residue_gemmi("//A/101", model);

        ASSERT_NE(res, nullptr);
        ASSERT_NE(upstream, nullptr);
        ASSERT_NE(downstream, nullptr);

        auto atoms = coot::residue_by_phi_psi::get_connecting_residue_atoms_gemmi(res, upstream, downstream);

        EXPECT_FALSE(atoms.empty());

        EXPECT_NEAR(atoms.N_pos.x(), 47.758, 1e-3);
        EXPECT_NEAR(atoms.N_pos.y(), 5.432, 1e-3);
        EXPECT_NEAR(atoms.N_pos.z(), 54.249, 1e-3);

        EXPECT_NEAR(atoms.CA_pos.x(), 48.663, 1e-3);
        EXPECT_NEAR(atoms.CA_pos.y(), 4.527, 1e-3);
        EXPECT_NEAR(atoms.CA_pos.z(), 53.538, 1e-3);

        EXPECT_NEAR(atoms.C_pos.x(), 47.897, 1e-3);
        EXPECT_NEAR(atoms.C_pos.y(), 3.657, 1e-3);
        EXPECT_NEAR(atoms.C_pos.z(), 52.543, 1e-3);

        EXPECT_TRUE(atoms.upstream_C.first);
        EXPECT_TRUE(atoms.downstream_N.first);

        EXPECT_NEAR(atoms.upstream_C.second.x(), 46.976, 1e-3);
        EXPECT_NEAR(atoms.upstream_C.second.y(), 6.304, 1e-3);
        EXPECT_NEAR(atoms.upstream_C.second.z(), 53.565, 1e-3);

        EXPECT_NEAR(atoms.downstream_N.second.x(), 48.493, 1e-3);
        EXPECT_NEAR(atoms.downstream_N.second.y(), 3.445, 1e-3);
        EXPECT_NEAR(atoms.downstream_N.second.z(), 51.359, 1e-3);

        auto phi = atoms.get_phi();
        auto psi = atoms.get_psi();

        EXPECT_TRUE(phi.first);
        EXPECT_TRUE(psi.first);
        EXPECT_NEAR(phi.second, -1.000638, 1e-4);
        EXPECT_NEAR(psi.second, 2.444335, 1e-4);
    }

    // Case 2: First residue — no upstream neighbour
    {
        gemmi::Residue *res = coot::get_residue_gemmi("//A/1", model);
        gemmi::Residue *downstream = coot::get_residue_gemmi("//A/2", model);

        ASSERT_NE(res, nullptr);
        ASSERT_NE(downstream, nullptr);

        auto atoms = coot::residue_by_phi_psi::get_connecting_residue_atoms_gemmi(res, nullptr, downstream);

        EXPECT_FALSE(atoms.empty());

        EXPECT_NEAR(atoms.N_pos.x(), 12.334, 1e-3);
        EXPECT_NEAR(atoms.N_pos.y(), 2.772, 1e-3);
        EXPECT_NEAR(atoms.N_pos.z(), 74.951, 1e-3);

        EXPECT_NEAR(atoms.CA_pos.x(), 13.401, 1e-3);
        EXPECT_NEAR(atoms.CA_pos.y(), 3.78, 1e-3);
        EXPECT_NEAR(atoms.CA_pos.z(), 74.703, 1e-3);

        EXPECT_NEAR(atoms.C_pos.x(), 13.976, 1e-3);
        EXPECT_NEAR(atoms.C_pos.y(), 4.253, 1e-3);
        EXPECT_NEAR(atoms.C_pos.z(), 76.035, 1e-3);

        EXPECT_FALSE(atoms.upstream_C.first);
        EXPECT_TRUE(atoms.downstream_N.first);
    }

    // Case 3: Last residue with no downstream
    {
        gemmi::Residue *res = coot::get_residue_gemmi("//A/224", model);
        gemmi::Residue *upstream = coot::get_residue_gemmi("//A/223", model);

        ASSERT_NE(res, nullptr);
        ASSERT_NE(upstream, nullptr);

        auto atoms = coot::residue_by_phi_psi::get_connecting_residue_atoms_gemmi(res, upstream, nullptr);

        EXPECT_FALSE(atoms.empty());

        EXPECT_NEAR(atoms.N_pos.x(), 25.949, 1e-3);
        EXPECT_NEAR(atoms.N_pos.y(), 9.076, 1e-3);
        EXPECT_NEAR(atoms.N_pos.z(), 23.107, 1e-3);

        EXPECT_NEAR(atoms.CA_pos.x(), 25.403, 1e-3);
        EXPECT_NEAR(atoms.CA_pos.y(), 7.776, 1e-3);
        EXPECT_NEAR(atoms.CA_pos.z(), 22.771, 1e-3);

        EXPECT_NEAR(atoms.C_pos.x(), 26.498, 1e-3);
        EXPECT_NEAR(atoms.C_pos.y(), 6.734, 1e-3);
        EXPECT_NEAR(atoms.C_pos.z(), 22.547, 1e-3);

        EXPECT_TRUE(atoms.upstream_C.first);
        EXPECT_FALSE(atoms.downstream_N.first);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}