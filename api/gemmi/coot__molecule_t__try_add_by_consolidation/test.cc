#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Test that consolidation works when source has single-residue chains
TEST(TryAddByConsolidation, SingleResidueChainConsolidation) {
    // Create target structure with one chain 'A' containing multiple residues
    gemmi::Structure target_st;
    target_st.name = "target";
    
    gemmi::Model target_model;
    gemmi::Chain target_chain("A");
    
    // Add residues to chain A (multiple residue types to prevent consolidation)
    gemmi::Residue res1;
    res1.name = "ALA";
    res1.seqid.num.value = 1;
    target_chain.residues.push_back(res1);
    
    gemmi::Residue res2;
    res2.name = "GLU";
    res2.seqid.num.value = 2;
    target_chain.residues.push_back(res2);
    
    target_model.chains.push_back(target_chain);
    target_st.models.push_back(target_model);
    
    // Create source structure with single-residue chain 'B' (should consolidate if match found)
    gemmi::Structure source_st;
    source_st.name = "source";
    
    gemmi::Model source_model;
    gemmi::Chain source_chain("B");
    
    // Chain B has only ALA residues - but we won't consolidate since target has mixed
    gemmi::Residue source_res1;
    source_res1.name = "ALA";
    source_res1.seqid.num.value = 10;
    source_chain.residues.push_back(source_res1);
    
    source_model.chains.push_back(source_chain);
    source_st.models.push_back(source_model);
    
    // Call the function under test
    std::pair<bool, std::vector<std::string>> result =
        coot::molecule_t_try_add_by_consolidation_gemmi(source_st, target_st);
    
    // The function should succeed (at least attempt consolidation)
    EXPECT_TRUE(result.first);
    
    // Should have processed at least one chain
    EXPECT_GE(result.second.size(), 1u);
}

// Test that chains with single residue type get consolidated when possible
TEST(TryAddByConsolidation, MultipleSingleResidueChains) {
    // Create target structure with chains that have single residue types
    gemmi::Structure target_st;
    target_st.name = "target";
    
    gemmi::Model target_model;
    
    // Chain 'A' has only CYS residues
    gemmi::Chain chain_a("A");
    for (int i = 1; i <= 10; ++i) {
        gemmi::Residue res;
        res.name = "CYS";
        res.seqid.num.value = i;
        chain_a.residues.push_back(res);
    }
    target_model.chains.push_back(chain_a);
    
    // Chain 'B' has only GLN residues
    gemmi::Chain chain_b("B");
    for (int i = 1; i <= 5; ++i) {
        gemmi::Residue res;
        res.name = "GLN";
        res.seqid.num.value = i;
        chain_b.residues.push_back(res);
    }
    target_model.chains.push_back(chain_b);
    
    target_st.models.push_back(target_model);
    
    // Create source structure with chains to consolidate
    gemmi::Structure source_st;
    source_st.name = "source";
    
    gemmi::Model source_model;
    
    // Chain 'X' has only CYS residues - should consolidate with chain A
    gemmi::Chain chain_x("X");
    for (int i = 20; i <= 25; ++i) {
        gemmi::Residue res;
        res.name = "CYS";
        res.seqid.num.value = i;
        chain_x.residues.push_back(res);
    }
    source_model.chains.push_back(chain_x);
    
    // Chain 'Y' has only GLN residues - should consolidate with chain B
    gemmi::Chain chain_y("Y");
    for (int i = 10; i <= 15; ++i) {
        gemmi::Residue res;
        res.name = "GLN";
        res.seqid.num.value = i;
        chain_y.residues.push_back(res);
    }
    source_model.chains.push_back(chain_y);
    
    source_st.models.push_back(source_model);
    
    // Call the function under test
    std::pair<bool, std::vector<std::string>> result =
        coot::molecule_t_try_add_by_consolidation_gemmi(source_st, target_st);
    
    // The function should succeed
    EXPECT_TRUE(result.first);
    
    // Should have consolidated 2 chains (X with A, Y with B)
    EXPECT_GE(result.second.size(), 2u);
    
    // Target should now have 2 chains (A and B) with consolidated residues
    EXPECT_EQ(target_st.models[0].chains.size(), 2u);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}