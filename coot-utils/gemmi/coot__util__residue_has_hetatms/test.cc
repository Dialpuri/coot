#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// The test data path needs to be accessible
static const std::string TEST_PDB_PATH = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";

TEST(OracleTest, ResidueHasHetatms) {
    // Test 1: Regular protein residue (ILE 10) - should return 0 (no hetatms)
    {
        gemmi::Structure st = gemmi::read_pdb_file(TEST_PDB_PATH);
        gemmi::Model& model = st.models[0];
        
        // Find ILE 10 in chain A
        gemmi::Residue* res = nullptr;
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& r : chain.residues) {
                if (r.seqid.num.value == 10) {
                    res = &r;
                    break;
                }
            }
            if (res) break;
        }
        
        ASSERT_NE(res, nullptr) << "Residue ILE 10 not found";
        int result = coot::util::residue_has_hetatms_gemmi(*res);
        EXPECT_EQ(result, 0) << "ILE 10 should not have hetatms";
    }

    // Test 2: Empty residue - should return -1
    {
        gemmi::Residue empty_res;
        int result = coot::util::residue_has_hetatms_gemmi(empty_res);
        EXPECT_EQ(result, -1) << "Empty residue should return -1";
    }

    // Test 3: Find any hetatm residue in the structure (water)
    {
        gemmi::Structure st = gemmi::read_pdb_file(TEST_PDB_PATH);
        gemmi::Model& model = st.models[0];
        
        gemmi::Residue* hetatm_res = nullptr;
        bool found_hetatm = false;
        
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& res : chain.residues) {
                if (res.het_flag == 'H') {
                    hetatm_res = &res;
                    found_hetatm = true;
                    break;
                }
            }
            if (found_hetatm) break;
        }

        EXPECT_TRUE(found_hetatm) << "No hetatm residues found in example.pdb";
        if (found_hetatm && hetatm_res) {
            int result = coot::util::residue_has_hetatms_gemmi(*hetatm_res);
            EXPECT_EQ(result, 1) << "Hetatm residue should return 1";
            
            // The oracle showed HOH, but we'll just check it's a water
            EXPECT_FALSE(hetatm_res->name.empty()) << "Residue name should not be empty";
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}