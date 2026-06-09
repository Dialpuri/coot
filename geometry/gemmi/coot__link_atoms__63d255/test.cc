#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, link_atoms_gemmi) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    ASSERT_FALSE(st.models.empty()) << "No models loaded";

    // Case 1: link with model_number = 1 (mimics model_p passed with GetSerNum()==1)
    {
        gemmi::Connection link;
        // Atom 1: chain A, res 1, no icode, atom CA, no altloc
        link.partner1.chain_name = "A";
        link.partner1.res_id.seqid = gemmi::SeqId(1, ' ');
        link.partner1.atom_name = "CA";
        link.partner1.altloc = '\0';

        // Atom 2: chain A, res 2, no icode, atom CA, no altloc
        link.partner2.chain_name = "A";
        link.partner2.res_id.seqid = gemmi::SeqId(2, ' ');
        link.partner2.atom_name = "CA";
        link.partner2.altloc = '\0';

        std::pair<coot::atom_spec_t, coot::atom_spec_t> result =
            coot::link_atoms_gemmi(link, 1);

        EXPECT_EQ(result.first.chain_id, "A");
        EXPECT_EQ(result.first.res_no, 1);
        EXPECT_EQ(result.first.atom_name, "CA");
        EXPECT_EQ(result.first.model_number, 1);
        EXPECT_EQ(result.second.chain_id, "A");
        EXPECT_EQ(result.second.res_no, 2);
        EXPECT_EQ(result.second.atom_name, "CA");
        EXPECT_EQ(result.second.model_number, 1);
    }

    // Case 2: link with model_number = 1 (mimics model_p = nullptr, atom_spec_t defaults to 1)
    {
        gemmi::Connection link;
        // Atom 1: chain B, res 50, no icode, atom N, no altloc
        link.partner1.chain_name = "B";
        link.partner1.res_id.seqid = gemmi::SeqId(50, ' ');
        link.partner1.atom_name = "N";
        link.partner1.altloc = '\0';

        // Atom 2: chain B, res 51, no icode, atom CA, no altloc
        link.partner2.chain_name = "B";
        link.partner2.res_id.seqid = gemmi::SeqId(51, ' ');
        link.partner2.atom_name = "CA";
        link.partner2.altloc = '\0';

        std::pair<coot::atom_spec_t, coot::atom_spec_t> result =
            coot::link_atoms_gemmi(link);  // default model_number = 1

        EXPECT_EQ(result.first.chain_id, "B");
        EXPECT_EQ(result.first.res_no, 50);
        EXPECT_EQ(result.first.atom_name, "N");
        EXPECT_EQ(result.first.model_number, 1);
        EXPECT_EQ(result.second.chain_id, "B");
        EXPECT_EQ(result.second.res_no, 51);
        EXPECT_EQ(result.second.atom_name, "CA");
        EXPECT_EQ(result.second.model_number, 1);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
