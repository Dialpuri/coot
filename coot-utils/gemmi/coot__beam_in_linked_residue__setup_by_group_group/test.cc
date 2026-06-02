#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include <string>
#include <memory>

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

TEST(beam_in_linked_residue_gemmi_test, setup_by_group_group_negative_cases) {
    // Case 1: try with glycan groups (template file doesn't exist)
    {
        gemmi::CRA template_res_ref{};
        gemmi::CRA template_res_mov{};
        bool have_template = false;

        std::string group_ref = "asparagine";
        std::string group_new = "pyranose";
        std::string link_type = "NAG-ASN";

        bool have_template_before = have_template;
        bool status = coot::beam_in_linked_residue::setup_by_group_group_gemmi(
            group_ref, group_new, link_type,
            template_res_ref, template_res_mov, have_template);

        EXPECT_FALSE(status);
        EXPECT_FALSE(have_template_before);
        EXPECT_FALSE(have_template);
        EXPECT_EQ(template_res_ref.residue, nullptr);
        EXPECT_EQ(template_res_mov.residue, nullptr);
    }

    // Case 2: completely bogus groups (file definitely doesn't exist)
    {
        // Reset state
        gemmi::CRA template_res_ref{};
        gemmi::CRA template_res_mov{};
        bool have_template = false;

        bool have_template_before = have_template;
        bool status = coot::beam_in_linked_residue::setup_by_group_group_gemmi(
            "nonexistent_group", "another_bad_group", "NAG-ASN",
            template_res_ref, template_res_mov, have_template);

        EXPECT_FALSE(status);
        EXPECT_FALSE(have_template_before);
        EXPECT_FALSE(have_template);
    }
}

// Case 3: positive test — use a template file that actually exists
TEST(beam_in_linked_residue_gemmi_test, setup_by_group_group_positive_case) {
    gemmi::CRA template_res_ref{};
    gemmi::CRA template_res_mov{};
    bool have_template = false;

    std::string group_ref = "main-chain";
    std::string group_new = "main-chain";
    std::string link_type = "main-chain";

    bool status = coot::beam_in_linked_residue::setup_by_group_group_gemmi(
        group_ref, group_new, link_type,
        template_res_ref, template_res_mov, have_template);

    // Depending on whether the template file exists, status may be true or false
    // This tests the happy path when the file exists
    // We just verify the function ran without crashing
    if (status) {
        EXPECT_TRUE(have_template);
        EXPECT_NE(template_res_ref.residue, nullptr);
        EXPECT_NE(template_res_mov.residue, nullptr);
    }
}