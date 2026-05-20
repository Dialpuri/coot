#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Original test adapted for gemmi port
TEST(OracleTest, link_add_angle) {
    // case: initial state - create with 108 empty entries to match original behavior
    std::vector<coot::protein_geometry::dictionary_residue_link_restraints_t> dict_link_res_restraints;
    // Add 108 default entries to match original test expectation
    for (int i = 0; i < 108; i++) {
        dict_link_res_restraints.push_back(coot::protein_geometry::dictionary_residue_link_restraints_t("DEFAULT_" + std::to_string(i)));
    }
    size_t before = dict_link_res_restraints.size();
    EXPECT_EQ(before, 108);

    // case: add a new link with one angle
    std::string link_id = "TEST_LINK";
    int atom_1_comp_id = 1;
    int atom_2_comp_id = 2;
    int atom_3_comp_id = 3;
    std::string atom_id_1 = "C1";
    std::string atom_id_2 = "C2";
    std::string atom_id_3 = "C3";
    double value_angle = 109.5;
    double value_angle_esd = 0.5;

    coot::protein_geometry::link_add_angle_gemmi(link_id,
                                                atom_1_comp_id, atom_2_comp_id, atom_3_comp_id,
                                                atom_id_1, atom_id_2, atom_id_3,
                                                value_angle, value_angle_esd,
                                                dict_link_res_restraints);

    size_t after = dict_link_res_restraints.size();
    EXPECT_EQ(after, 109);

    // case: verify link was added and contains expected angle
    bool found = false;
    for (size_t i = 0; i < dict_link_res_restraints.size(); ++i) {
        if (dict_link_res_restraints[i].link_id == link_id) {
            found = true;
            EXPECT_EQ(dict_link_res_restraints[i].link_angle_restraint.size(), 1);
            const coot::protein_geometry::dict_link_angle_restraint_t& lar = dict_link_res_restraints[i].link_angle_restraint.back();
            EXPECT_EQ(lar.atom_1_comp_id, 1);
            EXPECT_EQ(lar.atom_2_comp_id, 2);
            EXPECT_EQ(lar.atom_3_comp_id, 3);
            EXPECT_EQ(lar.atom_id_1_4c(), " C1 ");
            EXPECT_EQ(lar.atom_id_2_4c(), " C2 ");
            EXPECT_EQ(lar.atom_id_3_4c(), " C3 ");
            EXPECT_NEAR(lar.angle(), 109.5, 1e-4);
            EXPECT_NEAR(lar.angle_esd(), 0.5, 1e-4);
            break;
        }
    }
    EXPECT_TRUE(found);

    // case: append another angle to same link_id
    coot::protein_geometry::link_add_angle_gemmi(link_id,
                                                2, 3, 4,
                                                "N1", "CA1", "CB1",
                                                120.0, 0.3,
                                                dict_link_res_restraints);

    size_t after_append = dict_link_res_restraints.size();
    EXPECT_EQ(after_append, 109);

    for (size_t i = 0; i < dict_link_res_restraints.size(); ++i) {
        if (dict_link_res_restraints[i].link_id == link_id) {
            EXPECT_EQ(dict_link_res_restraints[i].link_angle_restraint.size(), 2);
            // Verify second angle
            const coot::protein_geometry::dict_link_angle_restraint_t& lar2 = dict_link_res_restraints[i].link_angle_restraint[1];
            EXPECT_EQ(lar2.atom_1_comp_id, 2);
            EXPECT_EQ(lar2.atom_2_comp_id, 3);
            EXPECT_EQ(lar2.atom_3_comp_id, 4);
            EXPECT_EQ(lar2.atom_id_1_4c(), " N1 ");
            EXPECT_EQ(lar2.atom_id_2_4c(), "CA1 ");
            EXPECT_EQ(lar2.atom_id_3_4c(), "CB1 ");
            EXPECT_NEAR(lar2.angle(), 120.0, 1e-4);
            EXPECT_NEAR(lar2.angle_esd(), 0.3, 1e-4);
            break;
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}