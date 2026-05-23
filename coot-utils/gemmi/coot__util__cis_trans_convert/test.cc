#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <vector>
#include "function.hh"

TEST(OracleTest, cis_trans_convert) {
    const char *pdb_path = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";
    
    // Load the test PDB using gemmi
    gemmi::Structure st = gemmi::read_pdb_file(pdb_path);
    
    // Get the first two residues from chain A
    gemmi::Model& model = st.models[0];
    ASSERT_FALSE(model.chains.empty());
    
    gemmi::Chain& chain = model.chains[0];
    ASSERT_GE(chain.residues.size(), 2);
    
    gemmi::Residue& res1 = chain.residues[0];
    gemmi::Residue& res2 = chain.residues[1];
    
    // Get initial atoms for logging - gemmi stores names without trailing spaces
    gemmi::Atom* at1_CA = nullptr;
    gemmi::Atom* at1_C = nullptr;
    gemmi::Atom* at2_CA = nullptr;
    gemmi::Atom* at2_N = nullptr;
    
    for (auto& atom : res1.atoms) {
        if (atom.name == "CA") at1_CA = &atom;
        if (atom.name == "C") at1_C = &atom;
    }
    
    for (auto& atom : res2.atoms) {
        if (atom.name == "CA") at2_CA = &atom;
        if (atom.name == "N") at2_N = &atom;
    }
    
    // Verify atoms exist
    ASSERT_NE(at1_CA, nullptr);
    ASSERT_NE(at1_C, nullptr);
    ASSERT_NE(at2_CA, nullptr);
    ASSERT_NE(at2_N, nullptr);

    // Store initial values
    float at1_CA_x = at1_CA->pos.x;
    float at1_CA_y = at1_CA->pos.y;
    float at1_CA_z = at1_CA->pos.z;
    float at1_C_x = at1_C->pos.x;
    float at1_C_y = at1_C->pos.y;
    float at1_C_z = at1_C->pos.z;
    float at2_CA_x = at2_CA->pos.x;
    float at2_CA_y = at2_CA->pos.y;
    float at2_CA_z = at2_CA->pos.z;
    float at2_N_x = at2_N->pos.x;
    float at2_N_y = at2_N->pos.y;
    float at2_N_z = at2_N->pos.z;
    
    // Load standard residues from the same PDB
    gemmi::Structure standard_st = gemmi::read_pdb_file(pdb_path);
    gemmi::Model& standard_model = standard_st.models[0];
    
    // Extract ALA residues as trans_residues
    std::vector<gemmi::Residue*> trans_residues;
    std::vector<gemmi::Residue*> cis_residues;
    
    for (auto& chain : standard_model.chains) {
        for (auto& residue : chain.residues) {
            if (residue.name == "ALA") {
                trans_residues.push_back(&residue);
            }
            if (residue.name == "ARG") {
                cis_residues.push_back(&residue);
            }
        }
    }
    
    ASSERT_GE(trans_residues.size(), 2);
    ASSERT_GE(cis_residues.size(), 2);
    
    // Prepare residue pair for conversion
    std::pair<gemmi::Residue*, gemmi::Residue*> mol_residues(&res1, &res2);
    
    // Call the converted function
    int istatus = coot::util::cis_trans_convert_gemmi(mol_residues, trans_residues, cis_residues);
    
    // After conversion
    // Expected values from oracle output
    EXPECT_EQ(istatus, 1);
    EXPECT_NEAR(at1_CA->pos.x, 16.6022, 1e-4);
    EXPECT_NEAR(at1_CA->pos.y, 5.22867, 1e-4);
    EXPECT_NEAR(at1_CA->pos.z, 72.5501, 1e-4);
    EXPECT_NEAR(at1_C->pos.x, 16.3662, 1e-4);
    EXPECT_NEAR(at1_C->pos.y, 6.60181, 1e-4);
    EXPECT_NEAR(at1_C->pos.z, 71.9261, 1e-4);
    EXPECT_NEAR(at2_CA->pos.x, -3.62645, 1e-4);
    EXPECT_NEAR(at2_CA->pos.y, 7.56676, 1e-4);
    EXPECT_NEAR(at2_CA->pos.z, 76.8122, 1e-4);
    EXPECT_NEAR(at2_N->pos.x, -2.48463, 1e-4);
    EXPECT_NEAR(at2_N->pos.y, 7.59393, 1e-4);
    EXPECT_NEAR(at2_N->pos.z, 77.7175, 1e-4);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}