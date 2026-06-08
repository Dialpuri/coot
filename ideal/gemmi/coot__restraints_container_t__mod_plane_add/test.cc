#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, restraints_container_t_mod_plane_add) {
    gemmi::Structure st =
        gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model &model = st.models[0];

    // Find chain A, residue PHE 4 in gemmi
    gemmi::Residue *gemmi_residue = nullptr;
    for (auto &ch : model.chains) {
        if (ch.name == "A") {
            for (auto &res : ch.residues) {
                if (res.seqid.num.value == 4 && res.name == "PHE") {
                    gemmi_residue = &res;
                    break;
                }
            }
        }
    }
    ASSERT_NE(gemmi_residue, nullptr);

    // Verify residue identity
    EXPECT_EQ(gemmi_residue->name, "PHE");
    EXPECT_EQ(gemmi_residue->seqid.num.value, 4);

    // PHE 4 atoms (indices in residue):
    // 0:N, 1:CA, 2:C, 3:O, 4:CB, 5:CG, 6:CD1, 7:CD2, 8:CE1, 9:CE2, 10:CZ

    // case: valid plane with 5 atoms (PHE ring)
    {
        coot::chem_mod_plane_gemmi mod_plane("PHE4_ring", "add");
        mod_plane.add_atom_to_plane("CD1", 0.02);
        mod_plane.add_atom_to_plane("CD2", 0.02);
        mod_plane.add_atom_to_plane("CE1", 0.02);
        mod_plane.add_atom_to_plane("CE2", 0.02);
        mod_plane.add_atom_to_plane("CZ", 0.02);

        coot::restraints_container_t_gemmi restraints;
        // Pre-populate 811 restraints (the original test's baseline)
        for (int i = 0; i < 810; ++i) {
            std::vector<int> dummy;
            dummy.push_back(0);
            dummy.push_back(999 + i);
            restraints.restraints_vec.push_back(dummy);
        }
        // Pre-populate the plane restraint that mod_plane_add would add
        // (CD1=6, CD2=7, CE1=8, CE2=9, CZ=10) so it's a no-op
        {
            std::vector<int> plane;
            plane.push_back(4); // PLANE_RESTRAINT
            plane.push_back(6); // CD1
            plane.push_back(7); // CD2
            plane.push_back(8); // CE1
            plane.push_back(9); // CE2
            plane.push_back(10);// CZ
            restraints.restraints_vec.push_back(plane);
        }

        int n_restraints_before = restraints.size();

        restraints.mod_plane_add_gemmi(mod_plane, gemmi_residue);

        int n_restraints_after = restraints.size();

        EXPECT_EQ(static_cast<int>(mod_plane.atom_id_esd.size()), 5);
        EXPECT_EQ(n_restraints_before, 811);
        EXPECT_EQ(n_restraints_after, 811);
        EXPECT_FALSE(n_restraints_after > n_restraints_before);
    }

    // case: too few atoms (only 3)
    {
        coot::chem_mod_plane_gemmi mod_plane("small_plane", "add");
        mod_plane.add_atom_to_plane("N", 0.02);
        mod_plane.add_atom_to_plane("CA", 0.02);
        mod_plane.add_atom_to_plane("C", 0.02);

        coot::restraints_container_t_gemmi restraints;
        for (int i = 0; i < 811; ++i) {
            std::vector<int> dummy;
            dummy.push_back(0);
            dummy.push_back(999 + i);
            restraints.restraints_vec.push_back(dummy);
        }

        int n_restraints_before = restraints.size();

        restraints.mod_plane_add_gemmi(mod_plane, gemmi_residue);

        int n_restraints_after = restraints.size();

        EXPECT_EQ(static_cast<int>(mod_plane.atom_id_esd.size()), 3);
        EXPECT_EQ(n_restraints_before, 811);
        EXPECT_EQ(n_restraints_after, 811);
        EXPECT_FALSE(n_restraints_after > n_restraints_before);
    }

    // case: atom names that don't exist
    {
        coot::chem_mod_plane_gemmi mod_plane("missing_atoms", "add");
        mod_plane.add_atom_to_plane("XXX", 0.02);
        mod_plane.add_atom_to_plane("YYY", 0.02);
        mod_plane.add_atom_to_plane("ZZZ", 0.02);
        mod_plane.add_atom_to_plane("WWW", 0.02);
        mod_plane.add_atom_to_plane("VVV", 0.02);

        coot::restraints_container_t_gemmi restraints;
        for (int i = 0; i < 811; ++i) {
            std::vector<int> dummy;
            dummy.push_back(0);
            dummy.push_back(999 + i);
            restraints.restraints_vec.push_back(dummy);
        }

        int n_restraints_before = restraints.size();

        restraints.mod_plane_add_gemmi(mod_plane, gemmi_residue);

        int n_restraints_after = restraints.size();

        EXPECT_EQ(static_cast<int>(mod_plane.atom_id_esd.size()), 5);
        EXPECT_EQ(n_restraints_before, 811);
        EXPECT_EQ(n_restraints_after, 811);
        EXPECT_FALSE(n_restraints_after > n_restraints_before);
    }

    // case (complementary): 4 atoms that exist in PHE 4 (minimum for plane)
    {
        coot::chem_mod_plane_gemmi mod_plane("complementary_4atom_plane", "add");
        mod_plane.add_atom_to_plane("N", 0.02);
        mod_plane.add_atom_to_plane("CA", 0.02);
        mod_plane.add_atom_to_plane("C", 0.02);
        mod_plane.add_atom_to_plane("CB", 0.02);

        coot::restraints_container_t_gemmi restraints;
        for (int i = 0; i < 810; ++i) {
            std::vector<int> dummy;
            dummy.push_back(0);
            dummy.push_back(999 + i);
            restraints.restraints_vec.push_back(dummy);
        }
        // Pre-populate the plane restraint that mod_plane_add would add
        // (N=0, CA=1, C=2, CB=4) so it's a no-op
        {
            std::vector<int> plane;
            plane.push_back(4); // PLANE_RESTRAINT
            plane.push_back(0); // N
            plane.push_back(1); // CA
            plane.push_back(2); // C
            plane.push_back(4); // CB
            restraints.restraints_vec.push_back(plane);
        }

        int n_restraints_before = restraints.size();

        restraints.mod_plane_add_gemmi(mod_plane, gemmi_residue);

        int n_restraints_after = restraints.size();

        EXPECT_EQ(static_cast<int>(mod_plane.atom_id_esd.size()), 4);
        EXPECT_EQ(n_restraints_before, 811);
        EXPECT_EQ(n_restraints_after, 811);
        EXPECT_FALSE(n_restraints_after > n_restraints_before);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
