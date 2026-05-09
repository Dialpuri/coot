#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, AddAlternativeConformation) {
    // Read the PDB file
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Get the first model
    gemmi::Model& model = st.models[0];
    
    // Test 1: Valid residue (ILE10 in chain A)
    {
        std::string cid = "//A/10";
        gemmi::Residue* res = coot::molecule_t::cid_to_residue(model, cid);
        
        ASSERT_NE(res, nullptr) << "Residue not found: " << cid;
        
        // Get atom count before
        int n_atoms_before = res->atoms.size();
        
        // Call the function
        int status = coot::molecule_t::add_alternative_conformation_gemmi(model, cid);
        
        // Get atom count after - need to find residue again since atoms were added
        gemmi::Residue* res_after = coot::molecule_t::cid_to_residue(model, cid);
        ASSERT_NE(res_after, nullptr);
        int n_atoms_after = res_after->atoms.size();
        
        // As per oracle output: status is 0 for success
        EXPECT_EQ(status, 0);
        // Atom count should double (ILE10: 8 -> 16)
        EXPECT_EQ(n_atoms_after, n_atoms_before * 2);
    }
    
    // Test 2: Invalid CID (non-existent residue)
    {
        std::string cid = "//X/999";
        int status = coot::molecule_t::add_alternative_conformation_gemmi(model, cid);
        // Function should return non-zero on failure
        EXPECT_NE(status, 0);
    }
    
    // Test 3: Valid residue in chain A (SER5 in chain A)
    {
        std::string cid = "//A/5";
        gemmi::Residue* res = coot::molecule_t::cid_to_residue(model, cid);
        
        ASSERT_NE(res, nullptr) << "Residue not found: " << cid;
        
        // Get atom count before
        int n_atoms_before = res->atoms.size();
        
        // Call the function
        int status = coot::molecule_t::add_alternative_conformation_gemmi(model, cid);
        
        // Get atom count after - need to find residue again since atoms were added
        gemmi::Residue* res_after = coot::molecule_t::cid_to_residue(model, cid);
        ASSERT_NE(res_after, nullptr);
        int n_atoms_after = res_after->atoms.size();
        
        // As per oracle output: status is 0 for success
        EXPECT_EQ(status, 0);
        // Atom count should double (SER5: 9 -> 18)
        EXPECT_EQ(n_atoms_after, n_atoms_before * 2);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}