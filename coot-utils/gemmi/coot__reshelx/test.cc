#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Test fixture path
static const std::string TEST_PDB_PATH = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";

TEST(ReshelxTest, ReshelxTest) {
    // case: normal PDB load
    {
        // Read PDB file into gemmi Structure
        gemmi::Structure st = gemmi::read_pdb_file(TEST_PDB_PATH);
        ASSERT_FALSE(st.models.empty()) << "load failed";
        
        // INPUT checks - get values before reshelx
        int chains_in = static_cast<int>(st.models[0].chains.size());
        size_t atoms_in = gemmi::count_atom_sites(st);
        
        // Get cell and spacegroup info
        bool has_cell_in = (st.cell.a > 0);
        std::string sg_in = st.spacegroup_hm;
        
        // Call reshelx_gemmi
        gemmi::Structure shelx_st = coot::reshelx_gemmi(st);
        
        // OUTPUT checks for Test 1
        EXPECT_NE(shelx_st.models.size(), 0u) << "shelx_st should be non-empty";
        
        if (!shelx_st.models.empty()) {
            int chains_out = static_cast<int>(shelx_st.models[0].chains.size());
            size_t atoms_out = gemmi::count_atom_sites(shelx_st);
            
            bool has_cell_out = (shelx_st.cell.a > 0);
            std::string sg_out = shelx_st.spacegroup_hm;
            
            // Assertions based on oracle output
            EXPECT_EQ(chains_out, 1);  // Output has 1 chain
            EXPECT_EQ(atoms_out, 2107);  // Output has 2107 atoms
            EXPECT_FALSE(has_cell_out);  // No cell in output
            EXPECT_NE(sg_out.size(), 0u);
            EXPECT_EQ(sg_out, "P 21 21 21");
            
            // Test 2: reshelx on already shelved structure
            gemmi::Structure shelx_st_2 = coot::reshelx_gemmi(shelx_st);
            EXPECT_NE(shelx_st_2.models.size(), 0u) << "reshelx on shelved structure should return non-empty";
            
            if (!shelx_st_2.models.empty()) {
                size_t atoms_double = gemmi::count_atom_sites(shelx_st_2);
                EXPECT_EQ(atoms_double, 2107);  // Should remain same
            }
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}