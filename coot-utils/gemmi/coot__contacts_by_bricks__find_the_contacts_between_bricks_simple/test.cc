#include <gtest/gtest.h>
#include <gemmi/mmread.hpp>
#include <gemmi/model.hpp>
#include <vector>
#include <set>
#include "function.hh"

TEST(OracleTest, contacts_by_bricks_find_contacts_simple) {
    // Load structure from CIF
    gemmi::Structure st = gemmi::read_structure_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");

    auto& model = st.models[0];

    // Collect all atom pointers and their residue IDs
    std::vector<const gemmi::Atom*> atom_ptrs;
    std::vector<unsigned int> atom_residue_ids;

    for (auto& chain : model.chains) {
        for (auto& residue : chain.residues) {
            unsigned int res_id = static_cast<unsigned int>(residue.seqid.num.value);
            for (auto& atom : residue.atoms) {
                atom_residue_ids.push_back(res_id);
                atom_ptrs.push_back(&atom);
            }
        }
    }

    ASSERT_FALSE(atom_ptrs.empty());
    EXPECT_EQ(static_cast<unsigned int>(atom_ptrs.size()), 2474u);

    std::set<unsigned int> fixed_indices;

    // Case 1: typical call with only_between_different_residues_flag=true
    {
        float dist_max = 4.0f;

        coot::contacts_by_bricks::contacts_by_bricks cbb =
            coot::contacts_by_bricks::contacts_by_bricks_gemmi(
                atom_ptrs.data(),
                static_cast<int>(atom_ptrs.size()),
                fixed_indices);
        cbb.dist_nbc_max = dist_max;

        int n_atoms_for_contacts = static_cast<int>(atom_ptrs.size());
        std::vector<std::set<unsigned int>> contacts(n_atoms_for_contacts);

        coot::contacts_by_bricks::find_the_contacts_between_bricks_simple_gemmi(
            cbb, &contacts, true, atom_residue_ids);

        // Count total contacts found
        unsigned int total_contacts = 0;
        for (size_t i = 0; i < contacts.size(); i++) {
            total_contacts += contacts[i].size();
        }
        EXPECT_EQ(total_contacts, 3468u);
    }

    // Case 2: include same-residue contacts
    {
        float dist_max = 4.0f;

        coot::contacts_by_bricks::contacts_by_bricks cbb =
            coot::contacts_by_bricks::contacts_by_bricks_gemmi(
                atom_ptrs.data(),
                static_cast<int>(atom_ptrs.size()),
                fixed_indices);
        cbb.dist_nbc_max = dist_max;

        int n_atoms_for_contacts = static_cast<int>(atom_ptrs.size());
        std::vector<std::set<unsigned int>> contacts(n_atoms_for_contacts);

        coot::contacts_by_bricks::find_the_contacts_between_bricks_simple_gemmi(
            cbb, &contacts, false, atom_residue_ids);

        unsigned int total_contacts = 0;
        for (size_t i = 0; i < contacts.size(); i++) {
            total_contacts += contacts[i].size();
        }
        EXPECT_EQ(total_contacts, 5502u);
    }

    // Case 3: very small distance (edge case - few/no contacts)
    {
        float dist_max = 0.5f;

        coot::contacts_by_bricks::contacts_by_bricks cbb =
            coot::contacts_by_bricks::contacts_by_bricks_gemmi(
                atom_ptrs.data(),
                static_cast<int>(atom_ptrs.size()),
                fixed_indices);
        cbb.dist_nbc_max = dist_max;

        int n_atoms_for_contacts = static_cast<int>(atom_ptrs.size());
        std::vector<std::set<unsigned int>> contacts(n_atoms_for_contacts);

        coot::contacts_by_bricks::find_the_contacts_between_bricks_simple_gemmi(
            cbb, &contacts, true, atom_residue_ids);

        unsigned int total_contacts = 0;
        for (size_t i = 0; i < contacts.size(); i++) {
            total_contacts += contacts[i].size();
        }
        EXPECT_EQ(total_contacts, 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}