#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, nmer_crankshaft_set) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    zo::rama_table_set zorts;

    // Case 1: n_peptides=3, mid residue A/10
    // Oracle output: cs_vec_size=0, residue_count=0
    {
        coot::residue_spec_t rs("A", 10);
        unsigned int n_peptides = 3;

        coot::nmer_crankshaft_set ncs(rs, n_peptides, zorts, st.models[0]);
        EXPECT_EQ(ncs.size(), 0u);

        std::vector<gemmi::CRA> res_vec = ncs.residues_gemmi();
        EXPECT_EQ(res_vec.size(), 0u);
    }

    // Case 2: empty nmer_crankshaft_set — edge case
    // Oracle output: empty_cs_vec_size=0, empty_residue_count=0
    {
        coot::nmer_crankshaft_set ncs_empty;
        EXPECT_EQ(ncs_empty.size(), 0u);
        std::vector<gemmi::CRA> res_vec = ncs_empty.residues_gemmi();
        EXPECT_EQ(res_vec.size(), 0u);
    }

    // Case 3: near end of chain A — res 290 with n_peptides=3
    // Oracle output: cs_vec_size=0, residue_count=0
    {
        coot::residue_spec_t rs("A", 290);
        unsigned int n_peptides = 3;

        coot::nmer_crankshaft_set ncs(rs, n_peptides, zorts, st.models[0]);
        EXPECT_EQ(ncs.size(), 0u);

        std::vector<gemmi::CRA> res_vec = ncs.residues_gemmi();
        EXPECT_EQ(res_vec.size(), 0u);
    }

    // Complementary case: explicit model_number=1 to pass get_residue guard
    // Chain A has residues 5-293 in model 1. Using a mid-chain residue
    // with explicit model=1. Observed behavior still produces 0 crankshafts
    // (likely due to rama_table_set not having entries for the residue types
    // or other internal conditions in make_one_crankshaft_set).
    {
        coot::residue_spec_t rs(1, "A", 50, "");
        unsigned int n_peptides = 3;

        EXPECT_NO_THROW({
            coot::nmer_crankshaft_set ncs(rs, n_peptides, zorts, st.models[0]);
            EXPECT_EQ(ncs.size(), 0u);

            std::vector<gemmi::CRA> res_vec = ncs.residues_gemmi();
            EXPECT_EQ(res_vec.size(), 0u);
        });
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}