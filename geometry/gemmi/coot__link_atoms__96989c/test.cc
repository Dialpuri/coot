#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/seqid.hpp>
#include "function.hh"

TEST(OracleTest, link_atoms) {
    // Load the test PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_EQ(st.models.size(), 1u);

    gemmi::Model* model = &st.models[0];

    // Check there are no connections (links) in the file
    EXPECT_EQ(st.connections.size(), 0u);

    // Case 3: synthetic Connection with model set
    {
        gemmi::Connection conn{};

        // partner1 fields
        conn.partner1.chain_name  = "2";
        conn.partner1.res_id.seqid.num = -2147483647;
        conn.partner1.res_id.seqid.icode = ' ';
        conn.partner1.atom_name = "ALA ";
        conn.partner1.altloc       = '\0';

        // partner2 fields
        conn.partner2.chain_name  = "5";
        conn.partner2.res_id.seqid.num = 31;
        conn.partner2.res_id.seqid.icode = ' ';
        conn.partner2.atom_name = "    ";
        conn.partner2.altloc       = '\0';

        auto pair = coot::link_atoms_gemmi(conn, model);
        EXPECT_EQ(pair.first.chain_id, "2");
        EXPECT_EQ(pair.first.res_no, -2147483647);
        EXPECT_EQ(pair.first.atom_name, "ALA ");
        EXPECT_EQ(pair.second.chain_id, "5");
        EXPECT_EQ(pair.second.res_no, 31);
        EXPECT_EQ(pair.second.atom_name, "    ");
        EXPECT_EQ(pair.first.model_number, 1);
    }

    // Case 4: synthetic Connection with nullptr model
    {
        gemmi::Connection conn{};

        conn.partner1.chain_name  = "2";
        conn.partner1.res_id.seqid.num = -2147483647;
        conn.partner1.res_id.seqid.icode = ' ';
        conn.partner1.atom_name = "ALA ";
        conn.partner1.altloc       = '\0';

        conn.partner2.chain_name  = "5";
        conn.partner2.res_id.seqid.num = 31;
        conn.partner2.res_id.seqid.icode = ' ';
        conn.partner2.atom_name = "    ";
        conn.partner2.altloc       = '\0';

        auto pair = coot::link_atoms_gemmi(conn, nullptr);
        EXPECT_EQ(pair.first.chain_id, "2");
        EXPECT_EQ(pair.first.res_no, -2147483647);
        EXPECT_EQ(pair.first.atom_name, "ALA ");
        EXPECT_EQ(pair.first.model_number, 1);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}