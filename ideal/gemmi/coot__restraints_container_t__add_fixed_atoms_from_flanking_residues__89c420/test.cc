#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, add_fixed_atoms_from_flanking_residues) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model& model = st.models[0];

    // Find chain A
    gemmi::Chain* chain = nullptr;
    for (gemmi::Chain& ch : model.chains) {
        if (ch.name == "A") {
            chain = &ch;
            break;
        }
    }
    ASSERT_NE(chain, nullptr);

    // Build restraints_container_t_gemmi
    // Collect non-H atoms from residues 1-20 (158 atoms) + residue 22 (11 atoms) = 169 total
    // This matches what the original make_restraints produces for this test fixture.
    coot::restraints_container_t_gemmi restraints;

    for (gemmi::Residue& res : chain->residues) {
        int seq = res.seqid.num.value;
        if (seq >= 1 && seq <= 20) {
            for (gemmi::Atom& atom : res.atoms) {
                if (!atom.element.is_hydrogen()) {
                    restraints.atom_residues.push_back(&res);
                }
            }
        } else if (seq == 22) {
            for (gemmi::Atom& atom : res.atoms) {
                if (!atom.element.is_hydrogen()) {
                    restraints.atom_residues.push_back(&res);
                }
            }
        }
    }

    restraints.n_atoms = restraints.atom_residues.size();

    // Record n_atoms (should be 169 as in original test)
    EXPECT_EQ(restraints.n_atoms, 169);

    // Case 1: Call with have_flanking_residue_at_start=true, have_flanking_residue_at_end=true
    {
        restraints.fixed_atom_indices.clear();

        size_t before = restraints.fixed_atom_indices.size();

        int iselection_start_res = 9;
        int iselection_end_res = 16;

        restraints.add_fixed_atoms_from_flanking_residues_gemmi(true, true,
                                                                 iselection_start_res,
                                                                 iselection_end_res);

        size_t after = restraints.fixed_atom_indices.size();

        EXPECT_EQ(before, 0U);
        EXPECT_EQ(after, 13U);
    }

    // Case 2: Call with both flags false — should not add anything
    {
        restraints.fixed_atom_indices.clear();

        size_t before = restraints.fixed_atom_indices.size();

        restraints.add_fixed_atoms_from_flanking_residues_gemmi(false, false, 9, 16);

        size_t after = restraints.fixed_atom_indices.size();

        EXPECT_EQ(before, 0U);
        EXPECT_EQ(after, 0U);
    }

    // Case 3: Call with only start flag true
    {
        restraints.fixed_atom_indices.clear();

        size_t before = restraints.fixed_atom_indices.size();

        restraints.add_fixed_atoms_from_flanking_residues_gemmi(true, false, 9, 16);

        size_t after = restraints.fixed_atom_indices.size();

        EXPECT_EQ(before, 0U);
        EXPECT_EQ(after, 9U);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}