#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, AtomTreeConstructor) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_GE(st.models.size(), 1) << "load failed";

    // Test 1: Create atom tree for a residue with valid contacts using constructor
    {
        gemmi::Chain& chain = st.models[0].chains[0];
        gemmi::Residue& res = chain.residues[0]; // first residue (0-indexed)

        EXPECT_EQ(res.seqid.num.value, 1) << "residue seqnum should be 1";

        // Get atoms from the residue
        const std::vector<gemmi::Atom>& atoms = res.atoms;
        int n_atoms = static_cast<int>(atoms.size());

        // First residue in example.pdb (ILE1) has 8 atoms
        EXPECT_EQ(n_atoms, 8) << "n_atoms_in_residue should be 8";

        if (n_atoms > 0) {
            // Create contact_indices with correct size (n_atoms)
            std::vector<std::vector<int>> contact_indices(n_atoms);

            // Create simple bonding pattern: atom i connects to atom i+1
            for (int i = 0; i < n_atoms - 1; i++) {
                contact_indices[i].push_back(i + 1);
                contact_indices[i + 1].push_back(i);
            }

            // Use constructor that should call fill_atom_vertex_vec_using_contacts_by_atom_selection
            coot::atom_tree_t tree;
            bool success = tree.fill_atom_vertex_vec_using_contacts_by_atom_selection_gemmi(contact_indices, atoms, 0);

            EXPECT_TRUE(success) << "fill_atom_vertex_vec_using_contacts_by_atom_selection should succeed";
            EXPECT_EQ(tree.atom_vertex_vec.size(), static_cast<size_t>(n_atoms)) << "atom_vertex_vec size should match n_atoms";
        }
    }

    // Test 2: Edge case - empty contact_indices (should fail early with return 0)
    {
        gemmi::Chain& chain = st.models[0].chains[0];
        gemmi::Residue& res = chain.residues[0]; // first residue (0-indexed)

        EXPECT_EQ(res.seqid.num.value, 1) << "residue seqnum should be 1 for edge case test";

        const std::vector<gemmi::Atom>& atoms = res.atoms;
        int n_atoms = static_cast<int>(atoms.size());

        std::vector<std::vector<int>> contact_indices; // empty - should trigger early return

        coot::atom_tree_t tree;
        bool success = tree.fill_atom_vertex_vec_using_contacts_by_atom_selection_gemmi(contact_indices, atoms, 0);

        EXPECT_FALSE(success) << "should fail early with return 0 for empty contact_indices";
        EXPECT_EQ(tree.atom_vertex_vec.size(), 0u) << "atom_vertex_vec size should be 0 for empty contact_indices";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}