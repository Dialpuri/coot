#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, ResidueRunAddGemmi) {
    // Load the PDB file using gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Case 1: Basic add with empty run
    {
        coot::util::map_to_model_correlation_stats_per_residue_run::residue_run_t run(5);
        EXPECT_EQ(run.residues.size(), 0u);
        
        // Collect some residues from the structure - residues 10, 11, 12
        std::vector<gemmi::Residue> rv;
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                for (auto& res : chain.residues) {
                    if (res.seqid.num.value == 10 || res.seqid.num.value == 11 || res.seqid.num.value == 12) {
                        rv.push_back(res);
                    }
                }
            }
        }
        
        run.add_gemmi(rv);
        
        EXPECT_EQ(run.residues.size(), 3u);
    }
    
    // Case 2: Add to non-empty run
    {
        coot::util::map_to_model_correlation_stats_per_residue_run::residue_run_t run(10);
        
        std::vector<gemmi::Residue> first;
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                for (auto& res : chain.residues) {
                    if (res.seqid.num.value == 2 || res.seqid.num.value == 3) {
                        first.push_back(res);
                    }
                }
            }
        }
        run.add_gemmi(first);
        EXPECT_EQ(run.residues.size(), 2u);
        
        std::vector<gemmi::Residue> second;
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                for (auto& res : chain.residues) {
                    if (res.seqid.num.value == 5 || res.seqid.num.value == 6 || res.seqid.num.value == 7) {
                        second.push_back(res);
                    }
                }
            }
        }
        run.add_gemmi(second);
        
        EXPECT_EQ(run.residues.size(), 5u);
    }
    
    // Case 3: Verify insert-at-beginning behavior
    {
        coot::util::map_to_model_correlation_stats_per_residue_run::residue_run_t run(10);
        
        std::vector<gemmi::Residue> first;
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                for (auto& res : chain.residues) {
                    if (res.seqid.num.value == 4 || res.seqid.num.value == 8) {
                        first.push_back(res);
                    }
                }
            }
        }
        run.add_gemmi(first);
        
        std::vector<gemmi::Residue> second;
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                for (auto& res : chain.residues) {
                    if (res.seqid.num.value == 9 || res.seqid.num.value == 10 || res.seqid.num.value == 11) {
                        second.push_back(res);
                    }
                }
            }
        }
        run.add_gemmi(second);
        
        EXPECT_EQ(run.residues.size(), 5u);
        EXPECT_EQ(run.residues[0].seqid.num.value, 9);
        EXPECT_EQ(run.residues[1].seqid.num.value, 10);
        EXPECT_EQ(run.residues[2].seqid.num.value, 11);
        EXPECT_EQ(run.residues[3].seqid.num.value, 4);
        EXPECT_EQ(run.residues[4].seqid.num.value, 8);
    }
    
    // Case 4: Add empty vector
    {
        coot::util::map_to_model_correlation_stats_per_residue_run::residue_run_t run(5);
        // Add first residue found in structure (residue 1)
        if (!st.models.empty() && !st.models[0].chains.empty() && !st.models[0].chains[0].residues.empty()) {
            run.residues.push_back(st.models[0].chains[0].residues[0]);
        }
        EXPECT_EQ(run.residues.size(), 1u);
        
        std::vector<gemmi::Residue> empty;
        run.add_gemmi(empty);
        
        EXPECT_EQ(run.residues.size(), 1u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}