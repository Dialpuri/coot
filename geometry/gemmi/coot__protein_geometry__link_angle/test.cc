#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/cifdoc.hpp>
#include "function.hh"

// Helper: creates a gemmi cif::Loop matching the original test data
static gemmi::cif::Loop make_test_loop() {
    gemmi::cif::Loop loop;
    loop.tags = {
        "link_id",
        "atom_id_1", "atom_id_2", "atom_id_3",
        "atom_1_comp_id", "atom_2_comp_id", "atom_3_comp_id",
        "value_angle", "value_angle_esd"
    };

    // Row 1: a valid angle entry
    loop.add_row({"test_link_1", "O", "C", "N", "1", "2", "3", "120.0", "2.0"});
    // Row 2: another valid angle entry
    loop.add_row({"test_link_2", "CA", "C", "O", "1", "2", "3", "110.5", "1.5"});

    return loop;
}

// Helper: fill a vector with 108 dummy entries (simulating init_standard())
static void fill_dummy_entries(std::vector<coot::protein_geometry::dictionary_residue_link_restraints_t>& v) {
    v.clear();
    for (int i = 0; i < 108; i++) {
        v.push_back(coot::protein_geometry::dictionary_residue_link_restraints_t("dummy"));
    }
}

TEST(OracleTest, link_angle) {
    {
        // Main case: from oracle — two valid angle entries added to loop
        std::vector<coot::protein_geometry::dictionary_residue_link_restraints_t> dict_link_res_restraints;
        fill_dummy_entries(dict_link_res_restraints);

        // Record state before calling link_angle
        size_t before = dict_link_res_restraints.size();

        // Build a gemmi cif::Loop for _chem_link_angle
        gemmi::cif::Loop loop = make_test_loop();

        EXPECT_EQ(loop.length(), 2u);
        EXPECT_EQ(loop.val(0, 0), "test_link_1");
        EXPECT_EQ(loop.val(1, 0), "test_link_2");

        // Call the function
        coot::protein_geometry::link_angle_gemmi(loop, dict_link_res_restraints);

        size_t after = dict_link_res_restraints.size();

        EXPECT_EQ(before, 108u);
        EXPECT_EQ(after, 110u);

        // Check the two newly added link_ids
        EXPECT_EQ(dict_link_res_restraints[108].link_id, "test_link_1");
        EXPECT_EQ(dict_link_res_restraints[109].link_id, "test_link_2");
    }

    // Complementary case: empty loop (exercises the early-return branch where
    // loop.GetLoopLength() == 0 and no rows are parsed)
    {
        std::vector<coot::protein_geometry::dictionary_residue_link_restraints_t> dict_link_res_restraints;
        fill_dummy_entries(dict_link_res_restraints);

        size_t before = dict_link_res_restraints.size();

        gemmi::cif::Loop empty_loop;
        empty_loop.tags = {"link_id"};
        // No rows added

        EXPECT_EQ(empty_loop.length(), 0u);

        EXPECT_NO_THROW(
            coot::protein_geometry::link_angle_gemmi(empty_loop, dict_link_res_restraints)
        );

        size_t after = dict_link_res_restraints.size();

        EXPECT_EQ(after, before);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}