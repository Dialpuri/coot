#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// ... TEST(...) blocks here ...

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

TEST(OracleTest, add_planar_peptide_restraint) {
    // Case: measure plane restraints for TRANS and PTRANS before and after add_planar_peptide_restraint
    {
        // Simulate init_standard() by creating the initial state that init_standard produces
        // init_standard() calls init(), which clears dict_link_res_restraints and planar_peptide_restraint_state
        // Then it reads default geometry files. From the test we know:
        //   - Before: TRANS has 2 plane restraints, PTRANS has 2 plane restraints
        //   - After: TRANS has 3, PTRANS has 3

        std::vector<coot::protein_geometry::dictionary_residue_link_restraints_t> dict_link_res_restraints;
        bool planar_peptide_restraint_state = false;

        // Simulate the initial state from init_standard() for TRANS link (2 plane restraints)
        {
            dict_link_res_restraints.push_back(
                coot::protein_geometry::dictionary_residue_link_restraints_t("TRANS"));
            // Add 2 initial plane restraints (from default geometry files)
            auto& trans = dict_link_res_restraints.back();
            trans.link_plane_restraint.push_back(
                coot::protein_geometry::dict_link_plane_restraint_t("plane-3-atoms", 0.08));
            trans.link_plane_restraint.back().atom_ids.push_back("CA");
            trans.link_plane_restraint.back().atom_comp_ids.push_back(1);
            trans.link_plane_restraint.back().atom_ids.push_back("C");
            trans.link_plane_restraint.back().atom_comp_ids.push_back(1);
            trans.link_plane_restraint.back().atom_ids.push_back("N");
            trans.link_plane_restraint.back().atom_comp_ids.push_back(2);

            trans.link_plane_restraint.push_back(
                coot::protein_geometry::dict_link_plane_restraint_t("plane-4-atoms", 0.08));
            trans.link_plane_restraint.back().atom_ids.push_back("C");
            trans.link_plane_restraint.back().atom_comp_ids.push_back(1);
            trans.link_plane_restraint.back().atom_ids.push_back("O");
            trans.link_plane_restraint.back().atom_comp_ids.push_back(1);
            trans.link_plane_restraint.back().atom_ids.push_back("N");
            trans.link_plane_restraint.back().atom_comp_ids.push_back(2);
            trans.link_plane_restraint.back().atom_ids.push_back("CA");
            trans.link_plane_restraint.back().atom_comp_ids.push_back(2);
        }

        // Simulate the initial state from init_standard() for PTRANS link (2 plane restraints)
        {
            dict_link_res_restraints.push_back(
                coot::protein_geometry::dictionary_residue_link_restraints_t("PTRANS"));
            auto& ptrans = dict_link_res_restraints.back();
            ptrans.link_plane_restraint.push_back(
                coot::protein_geometry::dict_link_plane_restraint_t("plane-3-atoms", 0.08));
            ptrans.link_plane_restraint.back().atom_ids.push_back("CA");
            ptrans.link_plane_restraint.back().atom_comp_ids.push_back(1);
            ptrans.link_plane_restraint.back().atom_ids.push_back("C");
            ptrans.link_plane_restraint.back().atom_comp_ids.push_back(1);
            ptrans.link_plane_restraint.back().atom_ids.push_back("N");
            ptrans.link_plane_restraint.back().atom_comp_ids.push_back(2);

            ptrans.link_plane_restraint.push_back(
                coot::protein_geometry::dict_link_plane_restraint_t("plane-4-atoms", 0.08));
            ptrans.link_plane_restraint.back().atom_ids.push_back("C");
            ptrans.link_plane_restraint.back().atom_comp_ids.push_back(1);
            ptrans.link_plane_restraint.back().atom_ids.push_back("O");
            ptrans.link_plane_restraint.back().atom_comp_ids.push_back(1);
            ptrans.link_plane_restraint.back().atom_ids.push_back("N");
            ptrans.link_plane_restraint.back().atom_comp_ids.push_back(2);
            ptrans.link_plane_restraint.back().atom_ids.push_back("CA");
            ptrans.link_plane_restraint.back().atom_comp_ids.push_back(2);
        }

        // Measure plane restraints for TRANS link BEFORE calling
        auto find_link = [&](const std::string& id) -> coot::protein_geometry::dictionary_residue_link_restraints_t* {
            for (auto& lr : dict_link_res_restraints) {
                if (lr.link_id == id) return &lr;
            }
            return nullptr;
        };

        auto trans_before = find_link("TRANS");
        int trans_plane_before = static_cast<int>(trans_before->link_plane_restraint.size());

        auto ptrans_before = find_link("PTRANS");
        int ptrans_plane_before = static_cast<int>(ptrans_before->link_plane_restraint.size());

        EXPECT_EQ(trans_plane_before, 2);
        EXPECT_EQ(ptrans_plane_before, 2);

        // Call the function
        coot::protein_geometry::add_planar_peptide_restraint_gemmi(dict_link_res_restraints, planar_peptide_restraint_state);

        // Measure plane restraints for TRANS and PTRANS link AFTER calling
        auto trans_after = find_link("TRANS");
        int trans_plane_after = static_cast<int>(trans_after->link_plane_restraint.size());

        auto ptrans_after = find_link("PTRANS");
        int ptrans_plane_after = static_cast<int>(ptrans_after->link_plane_restraint.size());

        EXPECT_EQ(trans_plane_after, 3);
        EXPECT_EQ(ptrans_plane_after, 3);
        EXPECT_EQ((trans_plane_after - trans_plane_before), 1);
        EXPECT_EQ((ptrans_plane_after - ptrans_plane_before), 1);
    }

    // Case: planar_peptide_restraint_state before and after add_planar_peptide_restraint
    {
        std::vector<coot::protein_geometry::dictionary_residue_link_restraints_t> dict_link_res_restraints;
        bool planar_peptide_restraint_state = false;

        bool state_before = planar_peptide_restraint_state;
        EXPECT_FALSE(state_before);

        coot::protein_geometry::add_planar_peptide_restraint_gemmi(dict_link_res_restraints, planar_peptide_restraint_state);
        bool state_after = planar_peptide_restraint_state;
        EXPECT_TRUE(state_after);
    }
}