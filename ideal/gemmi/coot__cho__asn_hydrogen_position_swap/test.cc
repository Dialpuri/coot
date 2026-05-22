#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <vector>
#include <utility>
#include <string>
#include "function.hh"

// Test the asn_hydrogen_position_swap_gemmi function
TEST(OracleTest, AsnHydrogenPositionSwap) {
    // Load test PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    ASSERT_FALSE(st.models.empty()) << "Failed to load test PDB";
    
    // Helper to find residue by chain/seq/icode
    auto find_res = [&st](const std::string& chain_id, int seqnum, char icode) -> gemmi::Residue* {
        for (gemmi::Model& model : st.models)
            for (gemmi::Chain& chain : model.chains)
                if (chain.name == chain_id)
                    for (gemmi::Residue& res : chain.residues)
                        if (res.seqid.num.value == seqnum && res.seqid.icode == icode)
                            return &res;
        return nullptr;
    };
    
    // Test case 1: Both residues null (guard clause test)
    {
        std::vector<std::pair<bool, gemmi::Residue*>> residues(2);
        residues[0] = std::make_pair(false, static_cast<gemmi::Residue*>(nullptr));
        residues[1] = std::make_pair(false, static_cast<gemmi::Residue*>(nullptr));
        
        // Function should not crash - early return works
        EXPECT_NO_THROW(coot::cho::asn_hydrogen_position_swap_gemmi(residues));
    }
    
    // Test case 2: First null, second valid (guard clause test)
    {
        gemmi::Residue* valid_res = find_res("A", 1, ' ');
        
        std::vector<std::pair<bool, gemmi::Residue*>> residues(2);
        residues[0] = std::make_pair(false, static_cast<gemmi::Residue*>(nullptr));
        residues[1] = std::make_pair(false, valid_res);
        
        // Function should not crash - early return works
        EXPECT_NO_THROW(coot::cho::asn_hydrogen_position_swap_gemmi(residues));
    }
    
    // Test case 3: Both valid (for ASN-NAG pair handling)
    {
        gemmi::Residue* res1 = find_res("A", 1, ' ');
        gemmi::Residue* res2 = find_res("A", 2, ' ');
        
        std::vector<std::pair<bool, gemmi::Residue*>> residues(2);
        residues[0] = std::make_pair(false, res1);
        residues[1] = std::make_pair(false, res2);
        
        // Function should not crash - early return works for non-ASN-NAG pair
        EXPECT_NO_THROW(coot::cho::asn_hydrogen_position_swap_gemmi(residues));
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}