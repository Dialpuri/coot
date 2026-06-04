#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, contact_info) {
    // Setup: load PDB
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        ASSERT_FALSE(st.models.empty());

        gemmi::Model& model = st.models[0];
        ASSERT_FALSE(model.chains.empty());

        gemmi::Chain& chain = model.chains[0];
        ASSERT_FALSE(chain.residues.empty());

        gemmi::Residue& res = chain.residues[0];
        ASSERT_GE(res.atoms.size(), 2u);

        // Build atom selection table
        std::vector<gemmi::Atom*> atom_table;
        for (gemmi::Atom& atom : res.atoms) {
            atom_table.push_back(&atom);
        }
        EXPECT_GE(static_cast<int>(atom_table.size()), 2);

        // Case: normal contacts (2 contacts)
        coot::mmdb_contact_gemmi contacts_raw[2];

        contacts_raw[0].id1 = 0;
        contacts_raw[0].id2 = 1;
        contacts_raw[0].dist = 0.0;

        contacts_raw[1].id1 = 1;
        contacts_raw[1].id2 = 0;
        contacts_raw[1].dist = 0.0;

        int nc = 2;

        coot::contact_info ci = coot::contact_info::contact_info_gemmi(atom_table, contacts_raw, nc);

        EXPECT_EQ(ci.n_contacts(), 2);
        EXPECT_EQ(ci.contacts[0].id1, 0);
        EXPECT_EQ(ci.contacts[0].id2, 1);
        EXPECT_EQ(ci.contacts[1].id1, 1);
        EXPECT_EQ(ci.contacts[1].id2, 0);

        // Case: empty contacts (nc=0)
        coot::mmdb_contact_gemmi* empty_contacts = new coot::mmdb_contact_gemmi[0];
        coot::contact_info ci_empty = coot::contact_info::contact_info_gemmi(atom_table, empty_contacts, 0);
        EXPECT_EQ(ci_empty.n_contacts(), 0);

        delete[] empty_contacts;
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}