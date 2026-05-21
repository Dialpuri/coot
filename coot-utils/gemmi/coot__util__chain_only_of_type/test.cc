#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, chain_only_of_type) {
    // Load the test structure
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty()) << "Failed to load PDB file";

    // Case 1: Look for a chain with only ALA residues
    // Chain A has many residue types, so this should return NULL
    {
        std::string residue_type = "ALA";
        const gemmi::Chain* chain = coot::util::chain_only_of_type_gemmi(st, residue_type);
        EXPECT_EQ(chain, nullptr) << "Chain not found expected for ALA (chain A has mixed residues)";
    }

    // Case 2: Look for a chain with only HOH residues
    // Chain B has only 1 residue which is HOH
    {
        std::string residue_type = "HOH";
        const gemmi::Chain* chain = coot::util::chain_only_of_type_gemmi(st, residue_type);
        EXPECT_NE(chain, nullptr) << "Chain should be found for HOH";
        if (chain != nullptr) {
            EXPECT_EQ(chain->name, "B") << "Expected chain B for HOH";
            EXPECT_EQ(chain->residues.size(), 1) << "Expected 1 residue in chain B";
        }
    }

    // Case 3: Look for a chain with non-existent residue type
    {
        std::string residue_type = "XXX";
        const gemmi::Chain* chain = coot::util::chain_only_of_type_gemmi(st, residue_type);
        EXPECT_EQ(chain, nullptr) << "Chain should not be found for non-existent residue type XXX";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}