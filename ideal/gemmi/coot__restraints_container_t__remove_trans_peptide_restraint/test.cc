#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(RemoveTransPeptideRestraintGemmiTest, removes_matching_restraint) {
    // Load the PDB to get real residues to work with
    gemmi::Structure st =
        gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];

    // Find chain A residues 10 and 11
    gemmi::Residue* res10 = nullptr;
    gemmi::Residue* res11 = nullptr;

    for (gemmi::Chain& chain : model.chains) {
        if (chain.name != "A") continue;
        for (gemmi::Residue& res : chain.residues) {
            if (res.seqid.num.value == 10) res10 = &res;
            if (res.seqid.num.value == 11) res11 = &res;
        }
    }

    ASSERT_NE(res10, nullptr);
    ASSERT_NE(res11, nullptr);

    // Build atom_residues mapping: index → residue*
    // Indices 0,1 → res10, indices 2,3 → res11
    std::vector<gemmi::Residue*> atom_residues = {res10, res10, res11, res11};

    // Create a TRANS_PEPTIDE_RESTRAINT that references atoms 0,1,2,3
    // → residues (res10, res10, res11, res11)
    simple_restraint_gemmi restraint;
    restraint.restraint_type = TRANS_PEPTIDE_RESTRAINT_VAL;
    restraint.atom_index_1 = 0;
    restraint.atom_index_2 = 1;
    restraint.atom_index_3 = 2;
    restraint.atom_index_4 = 3;
    restraint.is_active = true;

    std::vector<simple_restraint_gemmi> restraints = {restraint};

    // Call: should match and close the restraint
    remove_trans_peptide_restraint_gemmi(restraints, atom_residues, res10, res11);

    EXPECT_FALSE(restraints[0].is_active);
}

TEST(RemoveTransPeptideRestraintGemmiTest, does_not_match_reversed_order) {
    gemmi::Structure st =
        gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];

    gemmi::Residue* res10 = nullptr;
    gemmi::Residue* res11 = nullptr;

    for (gemmi::Chain& chain : model.chains) {
        if (chain.name != "A") continue;
        for (gemmi::Residue& res : chain.residues) {
            if (res.seqid.num.value == 10) res10 = &res;
            if (res.seqid.num.value == 11) res11 = &res;
        }
    }

    ASSERT_NE(res10, nullptr);
    ASSERT_NE(res11, nullptr);

    // Same setup as above
    std::vector<gemmi::Residue*> atom_residues = {res10, res10, res11, res11};

    simple_restraint_gemmi restraint;
    restraint.restraint_type = TRANS_PEPTIDE_RESTRAINT_VAL;
    restraint.atom_index_1 = 0;
    restraint.atom_index_2 = 1;
    restraint.atom_index_3 = 2;
    restraint.atom_index_4 = 3;
    restraint.is_active = true;

    std::vector<simple_restraint_gemmi> restraints = {restraint};

    // Reversed order: (res11, res10) should NOT match (res10, res10, res11, res11)
    remove_trans_peptide_restraint_gemmi(restraints, atom_residues, res11, res10);

    EXPECT_TRUE(restraints[0].is_active);
}

TEST(RemoveTransPeptideRestraintGemmiTest, null_pointers_no_crash) {
    gemmi::Structure st =
        gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];

    gemmi::Residue* res10 = nullptr;
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name != "A") continue;
        for (gemmi::Residue& res : chain.residues) {
            if (res.seqid.num.value == 10) res10 = &res;
        }
    }
    ASSERT_NE(res10, nullptr);

    std::vector<gemmi::Residue*> atom_residues = {res10, res10, res10, res10};

    simple_restraint_gemmi restraint;
    restraint.restraint_type = TRANS_PEPTIDE_RESTRAINT_VAL;
    restraint.atom_index_1 = 0;
    restraint.atom_index_2 = 1;
    restraint.atom_index_3 = 2;
    restraint.atom_index_4 = 3;
    restraint.is_active = true;

    std::vector<simple_restraint_gemmi> restraints = {restraint};

    // Null first, valid second — should not crash
    EXPECT_NO_THROW(
        remove_trans_peptide_restraint_gemmi(restraints, atom_residues, nullptr, res10));
    EXPECT_TRUE(restraints[0].is_active);

    // Valid first, null second — should not crash
    EXPECT_NO_THROW(
        remove_trans_peptide_restraint_gemmi(restraints, atom_residues, res10, nullptr));
    EXPECT_TRUE(restraints[0].is_active);

    // Both null — should not crash
    EXPECT_NO_THROW(
        remove_trans_peptide_restraint_gemmi(restraints, atom_residues, nullptr, nullptr));
    EXPECT_TRUE(restraints[0].is_active);
}

TEST(RemoveTransPeptideRestraintGemmiTest, skips_non_trans_peptide_restraints) {
    // Only TRANS_PEPTIDE_RESTRAINT entries are affected
    simple_restraint_gemmi bond;
    bond.restraint_type = 0; // NOT a trans peptide restraint
    bond.atom_index_1 = 0;
    bond.atom_index_2 = 1;
    bond.atom_index_3 = 2;
    bond.atom_index_4 = 3;
    bond.is_active = true;

    std::vector<simple_restraint_gemmi> restraints = {bond};
    std::vector<gemmi::Residue*> atom_residues;

    remove_trans_peptide_restraint_gemmi(restraints, atom_residues, nullptr, nullptr);
    EXPECT_TRUE(restraints[0].is_active);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}