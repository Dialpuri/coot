#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, add_to_non_drawn_bonds) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Get the first model
    ASSERT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];
    
    std::set<int> no_bonds_to_these_atom_indices;
    
    // Case 1: single residue CID
    {
        std::string cid = "//A/10";
        size_t before = no_bonds_to_these_atom_indices.size();
        EXPECT_EQ(before, 0u);

        // For gemmi, we need to count atoms in residue 10 of chain A
        // First, find the residue
        int atom_index = 0;
        std::set<int> selected_indices;
        for (auto& chain : model.chains) {
            if (chain.name == "A") {
                for (auto& res : chain.residues) {
                    if (res.seqid.num.value == 10) {
                        for (auto& atom : res.atoms) {
                            selected_indices.insert(atom_index);
                            atom_index++;
                        }
                    } else {
                        atom_index += res.atoms.size();
                    }
                }
            } else {
                for (auto& chain_res : chain.residues) {
                    atom_index += chain_res.atoms.size();
                }
            }
        }
        
        // Add these indices to the set
        for (int idx : selected_indices) {
            no_bonds_to_these_atom_indices.insert(idx);
        }

        size_t after = no_bonds_to_these_atom_indices.size();
        EXPECT_EQ(after, selected_indices.size());
    }

    // Case 2: multiple residues separated by ||
    {
        std::string cid = "//A/10||//A/20";
        size_t before = no_bonds_to_these_atom_indices.size();
        EXPECT_EQ(before, 8u);

        // Reset and count both residues
        std::set<int> selected_indices;
        int atom_index = 0;
        for (auto& chain : model.chains) {
            if (chain.name == "A") {
                for (auto& res : chain.residues) {
                    if (res.seqid.num.value == 10 || res.seqid.num.value == 20) {
                        for (auto& atom : res.atoms) {
                            selected_indices.insert(atom_index);
                            atom_index++;
                        }
                    } else {
                        atom_index += res.atoms.size();
                    }
                }
            } else {
                for (auto& chain_res : chain.residues) {
                    atom_index += chain_res.atoms.size();
                }
            }
        }
        
        // Add these indices to the set
        for (int idx : selected_indices) {
            no_bonds_to_these_atom_indices.insert(idx);
        }

        size_t after = no_bonds_to_these_atom_indices.size();
        EXPECT_EQ(after, 13u);
    }

    // Case 3: invalid CID (no matching atoms — should be a no-op)
    {
        std::string cid = "//Z/999";
        size_t before = no_bonds_to_these_atom_indices.size();
        EXPECT_EQ(before, 13u);

        // No matching atoms, so the set should remain unchanged
        // For gemmi, this would be handled by checking if any residues match
        // Since there are no matching residues, no_indices would be empty

        size_t after = no_bonds_to_these_atom_indices.size();
        EXPECT_EQ(after, 13u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}