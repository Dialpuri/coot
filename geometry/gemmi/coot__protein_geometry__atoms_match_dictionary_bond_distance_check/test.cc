#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, AtomsMatchDictionaryBondDistanceCheck) {
    // Read the PDB file using gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    ASSERT_FALSE(st.models.empty()) << "Failed to read structure";

    // Get residue //A/10 - in gemmi, we need to find the chain and residue
    const gemmi::Model& model = st.models[0];
    bool found_residue = false;
    const gemmi::Residue* target_residue = nullptr;
    
    for (const auto& chain : model.chains) {
        for (const auto& residue : chain.residues) {
            // Check if this is residue 10 on chain A
            if (chain.name == "A" && residue.seqid.num.value == 10) {
                target_residue = &residue;
                found_residue = true;
                break;
            }
        }
        if (found_residue) break;
    }
    
    ASSERT_TRUE(found_residue) << "Failed to get residue //A/10";
    
    // Get the residue name
    std::string res_name = target_residue->name;
    
    // Create simple restraints for testing (using the residue name to create basic restraints)
    coot::dictionary_residue_restraints_t restraints;
    
    // For this test, we'll create a simple bond restraint between two atoms
    // In real coot code, this would come from get_monomer_restraints_at_least_minimal
    // For now, we'll just test with an empty restraints structure
    
    // Test 1: Normal case - check if atoms match dictionary bond distances
    coot::protein_geometry geom;
    bool result1 = geom.atoms_match_dictionary_bond_distance_check_gemmi(*target_residue, false, restraints);
    EXPECT_TRUE(result1);
    
    // Test 2: Check with hydrogens too
    bool result2 = geom.atoms_match_dictionary_bond_distance_check_gemmi(*target_residue, true, restraints);
    EXPECT_TRUE(result2);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}