#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/metadata.hpp>
#include "function.hh"

TEST(OracleTestGemmi, remove_wrong_cis_peptides) {
    // Case 1: baseline — no added cis peptide
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        size_t n_before = st.cispeps.size();

        EXPECT_EQ(n_before, 1);

        coot::util::remove_wrong_cis_peptides_gemmi(st);

        size_t n_after = st.cispeps.size();
        EXPECT_EQ(n_after, 1);
    }

    // Case 2: Add a fake cis peptide that won't match real coordinates, then verify removal
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        // Add a fake cis peptide that clearly doesn't exist in coords
        gemmi::CisPep fake;
        fake.partner_c.chain_name = "A";
        fake.partner_c.res_id.seqid.num.value = 9999;
        fake.partner_c.res_id.seqid.icode = ' ';
        fake.partner_c.res_id.name = "GLY";
        fake.partner_c.atom_name = "";

        fake.partner_n.chain_name = "A";
        fake.partner_n.res_id.seqid.num.value = 9998;
        fake.partner_n.res_id.seqid.icode = ' ';
        fake.partner_n.res_id.name = "PRO";
        fake.partner_n.atom_name = "";

        st.cispeps.push_back(fake);

        size_t n_before = st.cispeps.size();
        EXPECT_EQ(n_before, 2);

        coot::util::remove_wrong_cis_peptides_gemmi(st);

        size_t n_after = st.cispeps.size();
        EXPECT_EQ(n_after, 1);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
