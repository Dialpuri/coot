#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/metadata.hpp>
#include <gemmi/seqid.hpp>
#include "function.hh"

TEST(OracleTest, atom_overlaps_container_is_linked) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    ASSERT_GT(st.models.size(), 0);
    gemmi::Model* model = &st.models[0];

    // Get first atom: first chain, first residue, first atom
    ASSERT_GT(model->chains.size(), 0);
    gemmi::Chain& chain_A = model->chains[0];
    ASSERT_GT(chain_A.residues.size(), 0);

    gemmi::Residue& res1 = chain_A.residues[0];
    ASSERT_GT(res1.atoms.size(), 0);
    gemmi::Atom& at1 = res1.atoms[0];

    // Second atom: second residue in same chain (if exists) or same residue
    gemmi::Residue* res2 = (chain_A.residues.size() > 1) ? &chain_A.residues[1] : &res1;
    ASSERT_GT(res2->atoms.size(), 0);
    gemmi::Atom& at2 = res2->atoms[0];

    gemmi::CRA cra1{&chain_A, &res1, &at1};
    gemmi::CRA cra2{&chain_A, res2, &at2};

    coot::atom_overlaps_container_t container;

    // Case 1: No link between atoms -> should return false
    {
        bool linked = container.is_linked_gemmi(cra1, cra2, model, &st);
        EXPECT_FALSE(linked);
    }

    // Case 2: Null atom -> should return false
    {
        gemmi::CRA cra_null{&chain_A, &res1, nullptr};
        bool linked = container.is_linked_gemmi(cra_null, cra2, model, &st);
        EXPECT_FALSE(linked);
    }

    // Case 3: Create a link and test -> should return true
    {
        gemmi::Connection conn;
        conn.type = gemmi::Connection::Type::Covale;

        // Set partner1
        conn.partner1.chain_name = chain_A.name;
        conn.partner1.res_id = gemmi::ResidueId(res1);
        conn.partner1.atom_name = at1.name;
        conn.partner1.altloc = at1.altloc;

        // Set partner2
        conn.partner2.chain_name = chain_A.name;
        conn.partner2.res_id = gemmi::ResidueId(*res2);
        conn.partner2.atom_name = at2.name;
        conn.partner2.altloc = at2.altloc;

        st.connections.push_back(conn);

        ASSERT_EQ(st.connections.size(), 1);

        bool linked = container.is_linked_gemmi(cra1, cra2, model, &st);
        EXPECT_TRUE(linked);

        // Test reverse order
        bool linked_rev = container.is_linked_gemmi(cra2, cra1, model, &st);
        EXPECT_TRUE(linked_rev);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}