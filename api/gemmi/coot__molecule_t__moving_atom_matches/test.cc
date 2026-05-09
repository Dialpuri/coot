#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, moving_atom_matches) {
    // Load a PDB file using gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Get first model
    gemmi::Model& model = st.models[0];
    
    // Build selection of all atoms as CRAs
    std::vector<gemmi::CRA> selection;
    for (auto& chain : model.chains)
        for (auto& res : chain.residues)
            for (auto& atom : res.atoms)
                selection.push_back({&chain, &res, &atom});
    
    // Ensure selection is populated
    ASSERT_GT(selection.size(), 0u);
    
    // Test 1: Check if first atom matches itself at index 0
    {
        bool result = coot::molecule_t::moving_atom_matches_gemmi(selection, 0, selection[0]);
        EXPECT_TRUE(result);
    }
    
    // Test 2: Check if second atom matches itself at index 1 (if available)
    if (selection.size() > 1) {
        bool result = coot::molecule_t::moving_atom_matches_gemmi(selection, 1, selection[1]);
        EXPECT_TRUE(result);
    }
    
    // Test 3: Test with out-of-bounds index (should return false)
    {
        bool result = coot::molecule_t::moving_atom_matches_gemmi(selection, 999999, selection[0]);
        EXPECT_FALSE(result);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}