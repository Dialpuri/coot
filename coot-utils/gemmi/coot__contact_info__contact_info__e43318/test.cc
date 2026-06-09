#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__molecule_t__get_residue/gemmi/function.hh"

#include "geometry/protein-geometry.hh"

TEST(OracleTest, contact_info) {
    coot::protein_geometry geom;
    geom.init_standard();

    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];

    // --- Case 1: Valid residue (first residue of chain A) with empty alt_conf ---
    {
        gemmi::Residue* res = coot::get_residue_gemmi("//A/1", model);
        ASSERT_NE(res, nullptr);

        std::vector<gemmi::Atom*> atom_selection;
        for (gemmi::Atom& a : res->atoms) {
            atom_selection.push_back(&a);
        }
        int n_atoms = static_cast<int>(atom_selection.size());

        std::string alt_conf = "";
        std::string monomer_type = res->name;

        EXPECT_EQ(monomer_type, "MET");
        EXPECT_EQ(n_atoms, 8);

        coot::contact_info::contact_info_result_gemmi ci = coot::contact_info::contact_info_gemmi(atom_selection, alt_conf, monomer_type, 0, &geom);

        EXPECT_EQ(ci.n_contacts(), 7);
    }

    // --- Case 2: First residue of chain B ---
    {
        gemmi::Residue* res = coot::get_residue_gemmi("//B/1", model);
        ASSERT_NE(res, nullptr);

        std::vector<gemmi::Atom*> atom_selection;
        for (gemmi::Atom& a : res->atoms) {
            atom_selection.push_back(&a);
        }
        int n_atoms = static_cast<int>(atom_selection.size());

        std::string alt_conf = "";
        std::string monomer_type = res->name;

        EXPECT_EQ(monomer_type, "HOH");
        EXPECT_EQ(n_atoms, 1);

        coot::contact_info::contact_info_result_gemmi ci = coot::contact_info::contact_info_gemmi(atom_selection, alt_conf, monomer_type, 0, &geom);

        EXPECT_EQ(ci.n_contacts(), 0);
    }

    // --- Case 3: Invalid alt_conf ("Z") should filter out all atoms, yielding expected contacts ---
    {
        gemmi::Residue* res = coot::get_residue_gemmi("//A/1", model);
        ASSERT_NE(res, nullptr);

        std::vector<gemmi::Atom*> atom_selection;
        for (gemmi::Atom& a : res->atoms) {
            atom_selection.push_back(&a);
        }
        int n_atoms = static_cast<int>(atom_selection.size());

        std::string alt_conf = "Z";
        std::string monomer_type = res->name;

        EXPECT_EQ(monomer_type, "MET");
        EXPECT_EQ(n_atoms, 8);

        coot::contact_info::contact_info_result_gemmi ci = coot::contact_info::contact_info_gemmi(atom_selection, alt_conf, monomer_type, 0, &geom);

        EXPECT_EQ(ci.n_contacts(), 7);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}