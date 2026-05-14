#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <string>
#include "coot/api/coot-molecule.hh"
#include "function.hh"

TEST(OracleTest, PiStackingContainer) {
    // Setup: load structure using gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Get residues using gemmi API
    std::vector<void*> residues;
    if (!st.models.empty()) {
        gemmi::Model& model = st.models[0];
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& res : chain.residues) {
                residues.push_back(static_cast<void*>(&res));
            }
        }
    }
    
    ASSERT_GT(residues.size(), 0u) << "No residues found";
    
    // Pick a reference ligand residue — look for non-protein residue
    void *res_ref = nullptr;
    for (size_t i = 0; i < residues.size(); i++) {
        gemmi::Residue *res = static_cast<gemmi::Residue*>(residues[i]);
        std::string res_name = res->name;
        if (res_name != "ALA" && res_name != "CYS" && res_name != "ASP" &&
            res_name != "GLU" && res_name != "PHE" && res_name != "GLY" &&
            res_name != "HIS" && res_name != "ILE" && res_name != "LYS" &&
            res_name != "LEU" && res_name != "MET" && res_name != "ASN" &&
            res_name != "PRO" && res_name != "GLN" && res_name != "ARG" &&
            res_name != "SER" && res_name != "THR" && res_name != "VAL" &&
            res_name != "TRP" && res_name != "TYR") {
            res_ref = residues[i];
            break;
        }
    }
    
    // Fallback to first residue if no ligand found
    if (!res_ref) {
        res_ref = residues[0];
    }
    
    // Create dictionary restraints for the ligand
    gemmi::Residue *res_ref_real = static_cast<gemmi::Residue*>(res_ref);
    std::string comp_id = res_ref_real->name;
    coot::dictionary_residue_restraints_t restraints(comp_id, 0);
    
    // Test case 1: valid ligand + environment
    {
        std::vector<std::vector<std::string> > aromatic_ring_list;
        pli::pi_stacking_container_t psc(restraints, residues, res_ref, aromatic_ring_list);
        EXPECT_EQ(psc.stackings.size(), 0u);
    }
    
    // Test case 2: empty residues — init should handle gracefully
    {
        std::vector<void*> empty_residues;
        std::vector<std::vector<std::string> > aromatic_ring_list;
        pli::pi_stacking_container_t psc2(restraints, empty_residues, res_ref, aromatic_ring_list);
        EXPECT_EQ(psc2.stackings.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}