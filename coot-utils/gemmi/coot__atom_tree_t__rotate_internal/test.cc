#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include <cmath>
#include <vector>

TEST(OracleTest, atom_tree_t_rotate_internal) {
    // Load PDB using gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Get residue A/10 (ILE) - has CA, CB sidechain
    gemmi::Model &model = st.models[0];
    
    // Find chain A
    gemmi::Chain *chain_ptr = nullptr;
    for (gemmi::Chain &chain : model.chains) {
        if (chain.name == "A") {
            chain_ptr = &chain;
            break;
        }
    }
    ASSERT_NE(chain_ptr, nullptr);
    
    // Find residue with seqid 10
    gemmi::Residue *res_ptr = nullptr;
    for (gemmi::Residue &res : chain_ptr->residues) {
        if (res.seqid.num.value == 10) {
            res_ptr = &res;
            break;
        }
    }
    ASSERT_NE(res_ptr, nullptr);
    
    const std::vector<gemmi::Atom> &residue_atoms = res_ptr->atoms;
    int n_residue_atoms = static_cast<int>(residue_atoms.size());
    
    ASSERT_EQ(n_residue_atoms, 8);
    
    // Build a simple linear contact_indices (chain: atom0-atom1-atom2-...)
    std::vector<std::vector<int>> contact(n_residue_atoms);
    for (int i = 0; i < n_residue_atoms; i++) {
        if (i > 0) contact[i].push_back(i - 1);
        if (i < n_residue_atoms - 1) contact[i].push_back(i + 1);
    }
    
    // Create a mutable copy of atoms so we can modify them
    std::vector<gemmi::Atom> mutable_atoms = residue_atoms;
    
    // We don't actually use the tree for anything in the original - the method just
    // takes the indices and modifies atoms directly via the atom table.
    // But to match the original, let's create an atom_tree_t (no-op since we don't use contact info)
    coot::atom_tree_t tree;
    
    // Case 1: Rotate last atom 90 degrees around first-two-atom axis
    {
        int moving_idx = n_residue_atoms - 1; // last atom (7)
        gemmi::Atom &moving_atom = mutable_atoms[moving_idx];
        double orig_x = moving_atom.pos.x;
        double orig_y = moving_atom.pos.y;
        double orig_z = moving_atom.pos.z;
        
        const gemmi::Atom &at0 = mutable_atoms[0];
        const gemmi::Atom &at1 = mutable_atoms[1];
        
        gemmi::Vec3 dir(at1.pos.x - at0.pos.x,
                        at1.pos.y - at0.pos.y,
                        at1.pos.z - at0.pos.z);
        gemmi::Vec3 base_pos(at1.pos.x, at1.pos.y, at1.pos.z);
        double angle = M_PI / 2.0; // 90 degrees
        
        std::vector<coot::map_index_t> moving;
        coot::map_index_t mi(moving_idx);
        moving.push_back(mi);
        
        EXPECT_NEAR(orig_x, 27.069, 1e-3);
        EXPECT_NEAR(orig_y, 8.313, 1e-3);
        EXPECT_NEAR(orig_z, 67.745, 1e-3);
        
        tree.rotate_internal_gemmi(mutable_atoms, moving, dir, base_pos, angle);
        
        EXPECT_NEAR(moving_atom.pos.x, 25.3494, 1e-3);
        EXPECT_NEAR(moving_atom.pos.y, 12.1696, 1e-3);
        EXPECT_NEAR(moving_atom.pos.z, 64.2354, 1e-3);
    }
    
    // Case 2: Zero angle rotation — no-op
    // Note: atoms are already rotated from Case 1, so "before" matches Case 1's "after"
    {
        int moving_idx = n_residue_atoms - 1;
        gemmi::Atom &moving_atom = mutable_atoms[moving_idx];
        double orig_x = moving_atom.pos.x;
        double orig_y = moving_atom.pos.y;
        double orig_z = moving_atom.pos.z;
        
        const gemmi::Atom &at0 = mutable_atoms[0];
        const gemmi::Atom &at1 = mutable_atoms[1];
        
        gemmi::Vec3 dir(at1.pos.x - at0.pos.x,
                        at1.pos.y - at0.pos.y,
                        at1.pos.z - at0.pos.z);
        gemmi::Vec3 base_pos(at1.pos.x, at1.pos.y, at1.pos.z);
        double angle = 0.0;
        
        std::vector<coot::map_index_t> moving;
        coot::map_index_t mi(moving_idx);
        moving.push_back(mi);
        
        EXPECT_NEAR(orig_x, 25.3494, 1e-3);
        EXPECT_NEAR(orig_y, 12.1696, 1e-3);
        EXPECT_NEAR(orig_z, 64.2354, 1e-3);
        
        tree.rotate_internal_gemmi(mutable_atoms, moving, dir, base_pos, angle);
        
        EXPECT_NEAR(moving_atom.pos.x, 25.3494, 1e-3);
        EXPECT_NEAR(moving_atom.pos.y, 12.1696, 1e-3);
        EXPECT_NEAR(moving_atom.pos.z, 64.2354, 1e-3);
    }
    
    // Case 3: Empty moving list — no-op
    {
        int moving_idx = n_residue_atoms - 1;
        gemmi::Atom &moving_atom = mutable_atoms[moving_idx];
        double orig_x = moving_atom.pos.x;
        
        std::vector<coot::map_index_t> moving; // empty
        
        const gemmi::Atom &at0 = mutable_atoms[0];
        const gemmi::Atom &at1 = mutable_atoms[1];
        gemmi::Vec3 dir(at1.pos.x - at0.pos.x,
                        at1.pos.y - at0.pos.y,
                        at1.pos.z - at0.pos.z);
        gemmi::Vec3 base_pos(at1.pos.x, at1.pos.y, at1.pos.z);
        double angle = M_PI / 4.0;
        
        EXPECT_NEAR(orig_x, 25.3494, 1e-3);
        
        tree.rotate_internal_gemmi(mutable_atoms, moving, dir, base_pos, angle);
        
        EXPECT_NEAR(moving_atom.pos.x, 25.3494, 1e-3);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}