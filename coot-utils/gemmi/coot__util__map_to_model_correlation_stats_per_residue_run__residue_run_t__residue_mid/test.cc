#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <vector>
#include <string>
#include "function.hh"

TEST(OracleTest, ResidueRunTest) {
    // case: setup
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        ASSERT_GE(st.models.size(), 1) << "PDB load failed";

        ASSERT_GE(st.models[0].chains.size(), 1);
        const gemmi::Chain& chain = st.models[0].chains[0];
        EXPECT_EQ(chain.residues.size(), 267);
    }

    // case: test case 1 - sufficient residues (5/5)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        ASSERT_GE(st.models.size(), 1);

        gemmi::Model& model = st.models[0];
        gemmi::Chain& chain = model.chains[0];

        std::vector<gemmi::Residue> test_residues;
        for (size_t i = 0; i < 5 && i < chain.residues.size(); i++) {
            test_residues.push_back(chain.residues[i]);
        }

        coot::util::map_to_model_correlation_stats_per_residue_run::residue_run_t run(5);
        for (const auto& res : test_residues) {
            run.add_residue(res);
        }

        const gemmi::Residue* mid_res = run.residue_mid();
        EXPECT_NE(mid_res, nullptr);
        
        if (mid_res) {
            EXPECT_EQ(mid_res->seqid.num.value, 3);
            EXPECT_EQ(mid_res->name, "ASN");
        }
    }

    // case: test case 2 - insufficient residues (3/5)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        ASSERT_GE(st.models.size(), 1);

        gemmi::Model& model = st.models[0];
        gemmi::Chain& chain = model.chains[0];

        std::vector<gemmi::Residue> test_residues;
        for (size_t i = 0; i < 3 && i < chain.residues.size(); i++) {
            test_residues.push_back(chain.residues[i]);
        }

        coot::util::map_to_model_correlation_stats_per_residue_run::residue_run_t run(5);
        for (const auto& res : test_residues) {
            run.add_residue(res);
        }

        const gemmi::Residue* mid_res = run.residue_mid();
        EXPECT_EQ(mid_res, nullptr);
    }

    // case: test case 3 - exact number (7/7)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        ASSERT_GE(st.models.size(), 1);

        gemmi::Model& model = st.models[0];
        gemmi::Chain& chain = model.chains[0];

        std::vector<gemmi::Residue> test_residues;
        for (size_t i = 0; i < 7 && i < chain.residues.size(); i++) {
            test_residues.push_back(chain.residues[i]);
        }

        coot::util::map_to_model_correlation_stats_per_residue_run::residue_run_t run(test_residues);
        const gemmi::Residue* mid_res = run.residue_mid();
        EXPECT_NE(mid_res, nullptr);
        
        if (mid_res) {
            EXPECT_EQ(mid_res->seqid.num.value, 4);
            EXPECT_EQ(mid_res->name, "PHE");
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}