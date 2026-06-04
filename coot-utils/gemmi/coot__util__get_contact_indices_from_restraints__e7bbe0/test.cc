#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "api/molecules-container.hh"
#include "geometry/protein-geometry.hh"

TEST(OracleTest, get_contact_indices_from_restraints) {
    molecules_container_t mc;
    mc.geometry_init_standard();

    int imol = mc.read_pdb("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_GT(imol, -1);

    coot::protein_geometry &geom = mc.get_geometry();

    // Use gemmi to load and get the residue
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());

    gemmi::Model &model = st.models[0];
    gemmi::Chain &chain = model.chains[0];
    gemmi::Residue &residue = chain.residues[0];

    ASSERT_NE(&residue, nullptr);
    EXPECT_EQ(residue.name, "MET");
    EXPECT_EQ(residue.seqid.num.value, 1);

    std::string resname = residue.name;
    int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;
    auto rp = geom.get_monomer_restraints(resname, imol_enc);
    ASSERT_TRUE(rp.first);

    coot::dictionary_residue_restraints_t &restraints = rp.second;

    EXPECT_EQ(residue.atoms.size(), 8u);

    // Case 1: regular_residue=true, add_reverse_contacts=false
    {
        bool regular = true;
        bool reverse = false;
        auto result = coot::util::get_contact_indices_from_restraints_gemmi(
            &residue, restraints, regular, reverse);

        EXPECT_EQ(result.size(), 8u);

        int total_contacts = 0;
        for (size_t i = 0; i < result.size(); i++) {
            total_contacts += result[i].size();
        }
        EXPECT_EQ(total_contacts, 7);

        // Verify first few contact entries
        EXPECT_EQ(result[0], std::vector<int>{1});
        EXPECT_EQ(result[1], (std::vector<int>{2, 4}));
        EXPECT_EQ(result[2], std::vector<int>{3});
        EXPECT_TRUE(result[3].empty());
        EXPECT_EQ(result[4], std::vector<int>{5});
    }

    // Case 2: regular_residue=true, add_reverse_contacts=true
    {
        bool regular = true;
        bool reverse = true;
        auto result = coot::util::get_contact_indices_from_restraints_gemmi(
            &residue, restraints, regular, reverse);

        EXPECT_EQ(result.size(), 8u);

        int total_contacts = 0;
        for (size_t i = 0; i < result.size(); i++) {
            total_contacts += result[i].size();
        }
        EXPECT_EQ(total_contacts, 14);
    }

    // Case 3: regular_residue=false (ligand mode), add_reverse_contacts=false
    {
        bool regular = false;
        bool reverse = false;
        auto result = coot::util::get_contact_indices_from_restraints_gemmi(
            &residue, restraints, regular, reverse);

        EXPECT_EQ(result.size(), 8u);

        int total_contacts = 0;
        for (size_t i = 0; i < result.size(); i++) {
            total_contacts += result[i].size();
        }
        EXPECT_EQ(total_contacts, 7);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}