#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(make_non_bonded_fixed_flags_gemmi_test, basic) {
    // Load PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());

    // Get chain A
    gemmi::Chain* chain = nullptr;
    for (gemmi::Model& model : st.models) {
        for (gemmi::Chain& ch : model.chains) {
            if (ch.name == "A") {
                chain = &ch;
                break;
            }
        }
    }
    ASSERT_NE(chain, nullptr);

    // Collect atoms and their residues from first few residues
    std::vector<gemmi::Atom> all_atoms;
    std::vector<gemmi::Residue*> residue_for_atom;

    for (gemmi::Residue& res : chain->residues) {
        for (gemmi::Atom& atom : res.atoms) {
            all_atoms.push_back(atom);
            residue_for_atom.push_back(&res);
        }
    }

    // Build restraints container
    coot::restraints_container_t_gemmi restraints;
    restraints.atom = all_atoms;
    restraints.residue_for_atom = residue_for_atom;
    restraints.n_atoms = all_atoms.size();

    // Use first 2 residues for controlled testing
    restraints.atom.clear();
    restraints.residue_for_atom.clear();
    restraints.fixed_atom_indices.clear();
    restraints.non_bonded_neighbour_residues.clear();

    auto res_it = chain->residues.begin();
    gemmi::Residue* res1 = &(*res_it);
    ++res_it;
    gemmi::Residue* res2 = &(*res_it);

    int idx = 0;
    for (gemmi::Atom& a : res1->atoms) {
        restraints.atom.push_back(a);
        restraints.residue_for_atom.push_back(res1);
        idx++;
    }
    for (gemmi::Atom& a : res2->atoms) {
        restraints.atom.push_back(a);
        restraints.residue_for_atom.push_back(res2);
        idx++;
    }
    restraints.n_atoms = idx;

    // Mark res2 as a non-bonded neighbour residue
    restraints.non_bonded_neighbour_residues.push_back(res2);

    // Mark some atoms from res1 as fixed
    if (res1->atoms.size() > 3) {
        restraints.fixed_atom_indices.insert(3);
    }

    // Case A: both atoms from res1 (not fixed, not in non_bonded residue)
    if (res1->atoms.size() >= 2) {
        std::vector<bool> result = restraints.make_non_bonded_fixed_flags_gemmi(0, 1);
        EXPECT_FALSE(result[0]);
        EXPECT_FALSE(result[1]);
    }

    // Case B: atom from res1 (not fixed) and atom from res2 (in non_bonded residue)
    int first_res2 = res1->atoms.size();
    if (res2->atoms.size() >= 1) {
        std::vector<bool> result = restraints.make_non_bonded_fixed_flags_gemmi(0, first_res2);
        EXPECT_FALSE(result[0]);  // res1 not in non_bonded
        EXPECT_TRUE(result[1]);   // res2 in non_bonded
    }

    // Case C: both indices 0 (same atom)
    std::vector<bool> result = restraints.make_non_bonded_fixed_flags_gemmi(0, 0);
    EXPECT_FALSE(result[0]);
    EXPECT_FALSE(result[1]);

    // Case D: a fixed atom and a non-bonded residue atom
    if (res1->atoms.size() > 3 && res2->atoms.size() >= 1) {
        std::vector<bool> result = restraints.make_non_bonded_fixed_flags_gemmi(3, first_res2);
        EXPECT_TRUE(result[0]);   // index 3 is fixed
        EXPECT_TRUE(result[1]);   // first_res2 in non_bonded residue
    }

    // Case E: both fixed atoms
    if (res1->atoms.size() > 5) {
        restraints.fixed_atom_indices.insert(5);
        std::vector<bool> result = restraints.make_non_bonded_fixed_flags_gemmi(3, 5);
        EXPECT_TRUE(result[0]);   // both in fixed
        EXPECT_TRUE(result[1]);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}