#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

static const std::string PDB_FILE = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";

TEST(OracleTest, merge_molecules_just_one_residue_at_given_spec) {
    // Load the main molecule
    gemmi::Structure mol = gemmi::read_pdb_file(PDB_FILE);
    coot::molecule_t mol_wrapper;
    mol_wrapper.structure = mol;
    
    // Test 1: Add valid single residue to new chain/position
    {
        // Find first residue in chain A
        gemmi::Residue* res = nullptr;
        for (auto& chain : mol_wrapper.model()->chains) {
            if (chain.name == "A" && !chain.residues.empty()) {
                res = &chain.residues.front();
                break;
            }
        }
        
        if (res) {
            std::string chain_id = "A";  // Fixed chain ID
            int seq_num = res->seqid.num.value;
            std::string ins_code(1, res->seqid.icode);
            
            // Create single residue molecule
            gemmi::Structure single_res_mol;
            gemmi::Model model;
            gemmi::Chain new_chain;
            new_chain.name = chain_id;
            gemmi::Residue res_copy = *res;
            
            // Copy atoms
            for (const auto& atom : res->atoms) {
                gemmi::Atom new_atom = atom;
                new_atom.pos = atom.pos;
                res_copy.atoms.push_back(new_atom);
            }
            
            new_chain.residues.push_back(res_copy);
            model.chains.push_back(new_chain);
            single_res_mol.models.push_back(model);
            
            coot::residue_spec_t target_spec("B", 999, "");
            
            int n_res_before = 0;
            for (const auto& m : mol_wrapper.model()->chains) {
                n_res_before += m.residues.size();
            }
            
            bool result = mol_wrapper.merge_molecules_just_one_residue_at_given_spec_gemmi(single_res_mol, target_spec);
            
            int n_res_after = 0;
            for (const auto& m : mol_wrapper.model()->chains) {
                n_res_after += m.residues.size();
            }
            
            EXPECT_TRUE(result);
            EXPECT_EQ(n_res_after, n_res_before + 1);
        }
    }
    
    // Test 2: Empty residue spec (should return false)
    {
        coot::residue_spec_t empty_spec;
        bool result2 = mol_wrapper.merge_molecules_just_one_residue_at_given_spec_gemmi(mol, empty_spec);
        EXPECT_FALSE(result2);
    }
    
    // Test 3: Target residue already exists (should return false)
    {
        coot::residue_spec_t existing_spec("A", 1, "");
        bool result3 = mol_wrapper.merge_molecules_just_one_residue_at_given_spec_gemmi(mol, existing_spec);
        EXPECT_FALSE(result3);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}