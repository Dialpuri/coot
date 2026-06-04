#include <gtest/gtest.h>
#include <gemmi/mmread.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include <set>
#include <vector>
#include <algorithm>

TEST(OracleTest, contacts_by_bricks) {
    // Load structure using gemmi
    gemmi::Structure st = gemmi::read_structure_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");
    ASSERT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];

    // Select all atoms
    std::vector<const gemmi::Atom*> valid_atoms;
    for (gemmi::Chain& chain : model.chains) {
        for (gemmi::Residue& residue : chain.residues) {
            for (gemmi::Atom& atom : residue.atoms) {
                valid_atoms.push_back(&atom);
            }
        }
    }

    int n_valid = static_cast<int>(valid_atoms.size());
    EXPECT_EQ(n_valid, 2474);
    const gemmi::Atom* atom_arr[] = {
        valid_atoms[0], valid_atoms[1], valid_atoms[2], valid_atoms[3],
        valid_atoms[4], valid_atoms[5], valid_atoms[6], valid_atoms[7],
        valid_atoms[8], valid_atoms[9],
    };

    // Case 1: all atoms, no fixed atoms
    {
        std::set<unsigned int> fixed_atom_indices;

        // Build an array of all atom pointers for gemmi
        std::vector<const gemmi::Atom*> atom_ptrs = valid_atoms;
        // We need to convert vector<const gemmi::Atom*> to const gemmi::Atom* const*
        // We can pass the data directly if we construct a proper array
        std::vector<const gemmi::Atom*> all_atom_ptrs = valid_atoms;

        auto contacts = coot::contacts_by_bricks::contacts_by_bricks_gemmi(
            all_atom_ptrs.data(), n_valid, fixed_atom_indices);

        EXPECT_EQ(contacts.n_atoms, 2474);
        EXPECT_EQ(contacts.brick_size, 20);
        EXPECT_EQ(contacts.dist_nbc_max, 8);
        EXPECT_EQ(contacts.atoms_in_bricks.size(), 48u);
        EXPECT_EQ(contacts.fixed_flags.size(), 2474u);
    }

    // Case 2: ligand LZA atoms (res 1299, chain A) with all marked as fixed
    {
        const gemmi::Residue* lza_res = nullptr;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue& res : chain.residues) {
                    if (res.seqid.num.value == 1299) {
                        lza_res = &res;
                        break;
                    }
                }
            }
        }

        std::vector<const gemmi::Atom*> lza_atoms;
        if (lza_res) {
            for (const gemmi::Atom& at : lza_res->atoms) {
                lza_atoms.push_back(&at);
            }
        }
        int n_lza = static_cast<int>(lza_atoms.size());
        EXPECT_EQ(n_lza, 25);

        const gemmi::Atom* lza_atom_arr = lza_atoms.empty() ? nullptr : lza_atoms.data()[0];

        std::set<unsigned int> fixed_atom_indices;
        for (int i = 0; i < n_lza; i++) {
            fixed_atom_indices.insert(static_cast<unsigned int>(i));
        }

        // Build array for the function
        std::vector<const gemmi::Atom*> lza_ptrs;
        for (int i = 0; i < n_lza; i++) {
            lza_ptrs.push_back(&(*lza_atoms[i]));
        }
        // Actually we need a const gemmi::Atom* const* — the lza_atoms vector gives us const gemmi::Atom* elements
        // lza_atoms.data() gives const gemmi::Atom**, which is what we need
        auto contacts = coot::contacts_by_bricks::contacts_by_bricks_gemmi(
            lza_atoms.data(), n_lza, fixed_atom_indices);

        EXPECT_EQ(contacts.n_atoms, 25);
        EXPECT_EQ(contacts.brick_size, 20);
        EXPECT_EQ(contacts.dist_nbc_max, 8);
        EXPECT_EQ(contacts.atoms_in_bricks.size(), 1u);
        EXPECT_EQ(contacts.fixed_flags.size(), 25u);
    }

    // Case 3: small subset - first 10 atoms, half fixed
    {
        int n_subset = std::min(n_valid, 10);
        std::vector<const gemmi::Atom*> subset_ptrs;
        for (int i = 0; i < n_subset; i++) {
            subset_ptrs.push_back(valid_atoms[i]);
        }

        std::set<unsigned int> fixed_atom_indices;
        fixed_atom_indices.insert(0);
        fixed_atom_indices.insert(2);
        fixed_atom_indices.insert(4);

        auto contacts = coot::contacts_by_bricks::contacts_by_bricks_gemmi(
            subset_ptrs.data(), n_subset, fixed_atom_indices);

        EXPECT_EQ(contacts.n_atoms, 10);
        EXPECT_EQ(contacts.brick_size, 20);
        EXPECT_EQ(contacts.dist_nbc_max, 8);
        EXPECT_EQ(contacts.atoms_in_bricks.size(), 1u);
        EXPECT_EQ(contacts.fixed_flags.size(), 10u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}