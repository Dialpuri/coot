#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/mmread.hpp>
#include "function.hh"

// Test for split_multi_model_molecule_gemmi function
TEST(OracleTest, split_multi_model_molecule) {
    // Case 1: Single-model structure (the provided test PDB)
    {
        std::string pdb_path = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif";
        gemmi::Structure st = gemmi::read_structure_file(pdb_path);
        
        EXPECT_EQ(st.models.size(), 1u) << "Input should have 1 model before splitting";
        
        std::vector<gemmi::Structure> result = coot::util::split_multi_model_molecule_gemmi(st);
        
        EXPECT_EQ(result.size(), 1u) << "Single-model structure should produce 1 result";
        
        // Verify each result has exactly 1 model
        for (const auto& split_st : result) {
            EXPECT_EQ(split_st.models.size(), 1u);
        }
    }
    
    // Case 2: Multi-model structure (create a simple one with 2 models)
    {
        gemmi::Structure st;
        
        // Add first model
        gemmi::Model model1(1);
        st.models.push_back(model1);
        
        // Add second model
        gemmi::Model model2(2);
        st.models.push_back(model2);
        
        EXPECT_EQ(st.models.size(), 2u) << "Input should have 2 models before splitting";
        
        std::vector<gemmi::Structure> result = coot::util::split_multi_model_molecule_gemmi(st);
        
        EXPECT_EQ(result.size(), 2u) << "Two-model structure should produce 2 results";
        
        // Verify each result has exactly 1 model
        bool all_single_model = true;
        for (size_t i = 0; i < result.size(); i++) {
            if (result[i].models.size() != 1u) {
                all_single_model = false;
            }
        }
        EXPECT_TRUE(all_single_model) << "Each split result should have exactly 1 model";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}