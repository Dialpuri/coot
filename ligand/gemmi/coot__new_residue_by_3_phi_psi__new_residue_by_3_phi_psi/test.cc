#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, new_residue_by_3_phi_psi) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];

    gemmi::Chain& chain = model.chains[0];
    EXPECT_LT(2, chain.residues.size());

    gemmi::Residue& res = chain.residues[10];
    gemmi::Residue& prev_res = chain.residues[9];

    {
        coot::new_residue_by_3_phi_psi nr3p("C", 
                                           {&chain, &res, nullptr},
                                           chain.residues);

        EXPECT_EQ(nr3p.chain_id, "A");
        EXPECT_EQ(nr3p.terminus_type, "C");
        EXPECT_EQ(nr3p.cra.residue->seqid.num.value, 11);
        EXPECT_EQ(nr3p.thread_pool_p, nullptr);
        EXPECT_EQ(nr3p.n_threads, 0);
        EXPECT_EQ(nr3p.upstream_neighbour_residue_p, nullptr);
        EXPECT_EQ(nr3p.downstream_neighbour_residue_p, nullptr);
    }

    {
        coot::new_residue_by_3_phi_psi nr3p("N",
                                           {&chain, &res, nullptr},
                                           chain.residues);

        EXPECT_EQ(nr3p.chain_id, "A");
        EXPECT_EQ(nr3p.terminus_type, "N");
        EXPECT_EQ(nr3p.cra.residue->seqid.num.value, 11);
        EXPECT_EQ(nr3p.thread_pool_p, nullptr);
        EXPECT_EQ(nr3p.n_threads, 0);
        EXPECT_EQ(nr3p.upstream_neighbour_residue_p, nullptr);
        EXPECT_EQ(nr3p.downstream_neighbour_residue_p, nullptr);
    }

    {
        coot::new_residue_by_3_phi_psi nr3p("C",
                                           {&chain, &res, nullptr},
                                           chain.residues);
        nr3p.upstream_neighbour_residue_p = &prev_res;

        EXPECT_EQ(nr3p.chain_id, "A");
        EXPECT_EQ(nr3p.terminus_type, "C");
        EXPECT_EQ(nr3p.cra.residue->seqid.num.value, 11);
        EXPECT_EQ(nr3p.thread_pool_p, nullptr);
        EXPECT_EQ(nr3p.n_threads, 0);
        EXPECT_NE(nr3p.upstream_neighbour_residue_p, nullptr);
        EXPECT_EQ(nr3p.upstream_neighbour_residue_p->name, "ILE");
        EXPECT_EQ(nr3p.upstream_neighbour_residue_p->seqid.num.value, 10);
        EXPECT_EQ(nr3p.downstream_neighbour_residue_p, nullptr);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}