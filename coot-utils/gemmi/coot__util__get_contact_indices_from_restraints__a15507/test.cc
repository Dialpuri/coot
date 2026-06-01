#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "api/molecules-container.hh"
#include "geometry/protein-geometry.hh"

TEST(OracleTest, get_contact_indices_from_restraints) {
    molecules_container_t mc;
    mc.geometry_init_standard();

    // Load PDB using gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    // Get protein geometry from molecules_container
    coot::protein_geometry *geom_p = &mc.get_geometry();

    // Get a residue: MET 1 from chain A (first model, first chain, first residue)
    gemmi::Model *model = &st.models[0];
    gemmi::Chain *chain = &model->chains[0];
    gemmi::Residue *residue = &chain->residues[0];

    ASSERT_NE(residue, nullptr);
    EXPECT_EQ(std::string(residue->name), "MET");
    EXPECT_EQ(residue->seqid.num.value, 1);
    EXPECT_EQ(residue->atoms.size(), 19u);

    // Case 1: Regular residue (MET 1), regular_residue_flag=true, add_reverse_contacts=false
    {
        bool regular_residue_flag = true;
        bool add_reverse_contacts = false;

        std::vector<std::vector<int>> contact_indices =
            coot::util::get_contact_indices_from_restraints_gemmi(residue, geom_p, regular_residue_flag, add_reverse_contacts);

        EXPECT_EQ(contact_indices.size(), 19u);

        EXPECT_EQ(contact_indices[0].size(), 4u);
        EXPECT_EQ(contact_indices[1].size(), 3u);
        EXPECT_EQ(contact_indices[2].size(), 1u);
        EXPECT_EQ(contact_indices[4].size(), 3u);
        EXPECT_EQ(contact_indices[5].size(), 3u);
        EXPECT_EQ(contact_indices[6].size(), 1u);
        EXPECT_EQ(contact_indices[7].size(), 3u);

        int total_contacts = 0;
        for (size_t i = 0; i < contact_indices.size(); i++) {
            total_contacts += contact_indices[i].size();
        }
        EXPECT_EQ(total_contacts, 18);
    }

    // Case 2: Regular residue (MET 1), regular_residue_flag=true, add_reverse_contacts=true
    {
        bool regular_residue_flag = true;
        bool add_reverse_contacts = true;

        std::vector<std::vector<int>> contact_indices =
            coot::util::get_contact_indices_from_restraints_gemmi(residue, geom_p, regular_residue_flag, add_reverse_contacts);

        EXPECT_EQ(contact_indices.size(), 19u);

        EXPECT_EQ(contact_indices[0].size(), 4u);
        EXPECT_EQ(contact_indices[1].size(), 4u);
        EXPECT_EQ(contact_indices[2].size(), 2u);
        EXPECT_EQ(contact_indices[3].size(), 1u);
        EXPECT_EQ(contact_indices[4].size(), 4u);
        EXPECT_EQ(contact_indices[5].size(), 4u);
        EXPECT_EQ(contact_indices[6].size(), 2u);
        EXPECT_EQ(contact_indices[7].size(), 4u);
        EXPECT_EQ(contact_indices[8].size(), 1u);
        EXPECT_EQ(contact_indices[9].size(), 1u);
        EXPECT_EQ(contact_indices[10].size(), 1u);
        EXPECT_EQ(contact_indices[11].size(), 1u);
        EXPECT_EQ(contact_indices[12].size(), 1u);
        EXPECT_EQ(contact_indices[13].size(), 1u);
        EXPECT_EQ(contact_indices[14].size(), 1u);
        EXPECT_EQ(contact_indices[15].size(), 1u);
        EXPECT_EQ(contact_indices[16].size(), 1u);
        EXPECT_EQ(contact_indices[17].size(), 1u);
        EXPECT_EQ(contact_indices[18].size(), 1u);

        int total_contacts = 0;
        for (size_t i = 0; i < contact_indices.size(); i++) {
            total_contacts += contact_indices[i].size();
        }
        EXPECT_EQ(total_contacts, 36);
    }

    // Case 3: Non-regular residue flag (MET 1), regular_residue_flag=false, add_reverse_contacts=false
    {
        bool regular_residue_flag = false;
        bool add_reverse_contacts = false;

        std::vector<std::vector<int>> contact_indices =
            coot::util::get_contact_indices_from_restraints_gemmi(residue, geom_p, regular_residue_flag, add_reverse_contacts);

        EXPECT_EQ(contact_indices.size(), 19u);

        EXPECT_EQ(contact_indices[0].size(), 0u);
        EXPECT_EQ(contact_indices[1].size(), 1u);
        EXPECT_EQ(contact_indices[2].size(), 1u);
        EXPECT_EQ(contact_indices[3].size(), 1u);
        EXPECT_EQ(contact_indices[4].size(), 1u);
        EXPECT_EQ(contact_indices[5].size(), 1u);
        EXPECT_EQ(contact_indices[6].size(), 1u);
        EXPECT_EQ(contact_indices[7].size(), 1u);
        EXPECT_EQ(contact_indices[8].size(), 1u);
        EXPECT_EQ(contact_indices[9].size(), 1u);
        EXPECT_EQ(contact_indices[10].size(), 1u);
        EXPECT_EQ(contact_indices[11].size(), 1u);
        EXPECT_EQ(contact_indices[12].size(), 1u);
        EXPECT_EQ(contact_indices[13].size(), 1u);
        EXPECT_EQ(contact_indices[14].size(), 1u);
        EXPECT_EQ(contact_indices[15].size(), 1u);
        EXPECT_EQ(contact_indices[16].size(), 1u);
        EXPECT_EQ(contact_indices[17].size(), 1u);
        EXPECT_EQ(contact_indices[18].size(), 1u);

        int total_contacts = 0;
        for (size_t i = 0; i < contact_indices.size(); i++) {
            total_contacts += contact_indices[i].size();
        }
        EXPECT_EQ(total_contacts, 18);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}