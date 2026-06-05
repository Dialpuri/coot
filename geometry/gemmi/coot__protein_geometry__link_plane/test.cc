#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/cifdoc.hpp>
#include "function.hh"

TEST(OracleTest, link_plane) {
    {
        // Case 1: First link_plane call
        // Simulate what coot::protein_geometry::init_standard() does for link restraints
        auto dict_link_res_restraints = coot::protein_geometry::init_standard_link_restraints();

        // Count link restraints BEFORE
        std::size_t before = dict_link_res_restraints.size();
        EXPECT_EQ(before, 108u);

        // Create an mmCIF loop with _chem_link_plane columns
        gemmi::cif::Loop loop;
        loop.tags = { "link_id", "atom_id", "atom_comp_id", "plane_id", "dist_esd" };
        loop.add_row({ "TEST_LINK", "A/CA/10", "1", "PLANE1", "0.5" });

        // Call link_plane
        coot::protein_geometry::link_plane_gemmi(loop, dict_link_res_restraints);

        // Count link restraints AFTER
        std::size_t after = dict_link_res_restraints.size();
        EXPECT_EQ(after, 109u);

        // Find and assert the TEST_LINK entry
        bool found = false;
        for (auto const &entry : dict_link_res_restraints) {
            if (entry.link_id == "TEST_LINK" && entry.link_plane_restraint.size() > 0) {
                found = true;
                auto const &pr = entry.link_plane_restraint[0];
                EXPECT_EQ(entry.link_id, "TEST_LINK");
                EXPECT_EQ(pr.plane_id, "PLANE1");
                EXPECT_EQ(pr.n_atoms(), 1u);
                EXPECT_EQ(pr.atom_id(0), "A/CA/10");
                EXPECT_EQ(pr.atom_comp_ids[0], 1);
                EXPECT_NEAR(pr.dist_esd(), 0.5, 1e-4);
                break;
            }
        }
        EXPECT_TRUE(found);

        // Case 2: Second link_plane call with different link_id
        gemmi::cif::Loop loop2;
        loop2.tags = { "link_id", "atom_id", "atom_comp_id", "plane_id", "dist_esd" };
        loop2.add_row({ "SECOND_LINK", "B/N/20", "2", "PLANE2", "1.0" });

        coot::protein_geometry::link_plane_gemmi(loop2, dict_link_res_restraints);

        std::size_t after2 = dict_link_res_restraints.size();
        EXPECT_EQ(after2, 110u);

        // Find and assert the SECOND_LINK entry
        bool found2 = false;
        for (auto const &entry : dict_link_res_restraints) {
            if (entry.link_id == "SECOND_LINK" && entry.link_plane_restraint.size() > 0) {
                found2 = true;
                auto const &pr = entry.link_plane_restraint[0];
                EXPECT_EQ(entry.link_id, "SECOND_LINK");
                EXPECT_EQ(pr.plane_id, "PLANE2");
                EXPECT_EQ(pr.n_atoms(), 1u);
                EXPECT_EQ(pr.atom_id(0), "B/N/20");
                EXPECT_EQ(pr.atom_comp_ids[0], 2);
                EXPECT_NEAR(pr.dist_esd(), 1.0, 1e-4);
                break;
            }
        }
        EXPECT_TRUE(found2);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}