#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, water_coordination_t) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty()) << "No models in PDB";

    // Case 1: typical radius (4.0 Angstrom) - should find water and contacts
    {
        double radius = 4.0;

        coot::util::water_coordination_t wc(st, 0, radius);

        std::vector<coot::util::contact_atoms_info_t> contacts = wc.get_contacts();
        EXPECT_EQ(contacts.size(), 1u);

        // Show details of first contact if present
        if (!contacts.empty()) {
            coot::util::contact_atoms_info_t& ci = contacts[0];
            gemmi::Atom const* central = ci.central_atom();
            ASSERT_NE(central, nullptr);
            EXPECT_EQ(ci.central_atom_name(), " O  ");
            ASSERT_NE(ci.cra.residue, nullptr);
            EXPECT_EQ(ci.cra.residue->name, "HOH");
            ASSERT_NE(ci.cra.chain, nullptr);
            EXPECT_EQ(ci.cra.chain->name, "B");
            EXPECT_EQ(ci.cra.residue->seqid.num.value, 1);
            EXPECT_EQ(ci.size(), 9u);
        }
    }

    // Case 2: very small radius (0.1 Angstrom) - should find water but no contacts
    {
        double radius = 0.1;

        coot::util::water_coordination_t wc2(st, 0, radius);

        std::vector<coot::util::contact_atoms_info_t> contacts2 = wc2.get_contacts();
        EXPECT_EQ(contacts2.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}