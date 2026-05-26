#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, SphericallyAveragedMolecule) {
    // Test case 1: Empty structure (no models)
    {
        gemmi::Structure st;
        float bin_size = 1.0;
        std::vector<std::pair<double, double>> result = 
            coot::util::spherically_averaged_molecule_gemmi(st, bin_size);
        
        EXPECT_EQ(result.size(), 0u);
    }
    
    // Test case 2: With actual PDB file
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        
        if (!st.models.empty()) {
            float bin_size = 2.0;
            std::vector<std::pair<double, double>> result = 
                coot::util::spherically_averaged_molecule_gemmi(st, bin_size);
            
            // Count total atoms
            size_t total_atoms = 0;
            for (const auto& model : st.models) {
                for (const auto& chain : model.chains) {
                    for (const auto& res : chain.residues) {
                        total_atoms += res.atoms.size();
                    }
                }
            }
            
            EXPECT_EQ(total_atoms, 2107u);
            EXPECT_EQ(result.size(), 23u);
            
            // Check first 5 bins
            EXPECT_EQ(result[0].first, 1);
            EXPECT_EQ(result[0].second, 1);
            EXPECT_EQ(result[1].first, 3);
            EXPECT_EQ(result[1].second, 9);
            EXPECT_EQ(result[2].first, 5);
            EXPECT_EQ(result[2].second, 31);
            EXPECT_EQ(result[3].first, 7);
            EXPECT_EQ(result[3].second, 73);
            EXPECT_EQ(result[4].first, 9);
            EXPECT_EQ(result[4].second, 87);
        } else {
            FAIL() << "no models in structure";
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}