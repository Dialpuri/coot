#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, LinkAddBond) {
    std::vector<coot::dictionary_residue_link_restraints_t> dict_link_res_restraints;
    for (int i = 0; i < 108; i++) {
        dict_link_res_restraints.push_back(coot::dictionary_residue_link_restraints_t("link_" + std::to_string(i)));
    }

    {
        size_t before_size = dict_link_res_restraints.size();
        EXPECT_EQ(before_size, 108);

        std::string link_id = "test_link";
        int atom_1_comp_id = 1;
        int atom_2_comp_id = 1;
        std::string atom_id_1 = "CA";
        std::string atom_id_2 = "CB";
        double value_dist = 1.5;
        double value_dist_esd = 0.02;

        coot::link_add_bond_gemmi(link_id, atom_1_comp_id, atom_2_comp_id,
                                  atom_id_1, atom_id_2, value_dist, value_dist_esd, dict_link_res_restraints);

        size_t after_size = dict_link_res_restraints.size();
        EXPECT_EQ(after_size, 109);

        bool found = false;
        size_t bond_count = 0;
        for (unsigned int i = 0; i < dict_link_res_restraints.size(); i++) {
            if (dict_link_res_restraints[i].link_id == link_id) {
                found = true;
                bond_count = dict_link_res_restraints[i].link_bond_restraint.size();
            }
        }
        EXPECT_TRUE(found);
        EXPECT_EQ(bond_count, 1);
    }

    {
        size_t before_size = dict_link_res_restraints.size();
        EXPECT_EQ(before_size, 109);

        std::string link_id = "test_link";
        int atom_1_comp_id = 1;
        int atom_2_comp_id = 1;
        std::string atom_id_1 = "N";
        std::string atom_id_2 = "C";
        double value_dist = 1.3;
        double value_dist_esd = 0.015;

        coot::link_add_bond_gemmi(link_id, atom_1_comp_id, atom_2_comp_id,
                                  atom_id_1, atom_id_2, value_dist, value_dist_esd, dict_link_res_restraints);

        size_t after_size = dict_link_res_restraints.size();
        EXPECT_EQ(after_size, 109);

        bool found = false;
        size_t bond_count = 0;
        for (unsigned int i = 0; i < dict_link_res_restraints.size(); i++) {
            if (dict_link_res_restraints[i].link_id == link_id) {
                found = true;
                bond_count = dict_link_res_restraints[i].link_bond_restraint.size();
            }
        }
        EXPECT_TRUE(found);
        EXPECT_EQ(bond_count, 2);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}