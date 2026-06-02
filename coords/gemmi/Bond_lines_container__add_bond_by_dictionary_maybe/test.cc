#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: find an atom and its residue by chain name, residue sequence number, and atom name
static const gemmi::Atom* find_atom_with_residue(
    const gemmi::Model& model,
    const std::string& chain_name,
    int seqnum,
    const std::string& atom_name,
    const gemmi::Residue** out_residue = nullptr)
{
    for (const auto& chain : model.chains) {
        if (chain.name != chain_name) continue;
        for (const auto& res : chain.residues) {
            if (res.seqid.num.value != seqnum) continue;
            for (const auto& atom : res.atoms) {
                if (atom.name == atom_name) {
                    if (out_residue) *out_residue = &res;
                    return &atom;
                }
            }
        }
    }
    return nullptr;
}

static gemmi::Residue* find_residue(
    gemmi::Model& model,
    const std::string& chain_name,
    int seqnum)
{
    for (auto& chain : model.chains) {
        if (chain.name != chain_name) continue;
        for (auto& res : chain.residues) {
            if (res.seqid.num.value == seqnum)
                return &res;
        }
    }
    return nullptr;
}

TEST(OracleTest, add_bond_by_dictionary_maybe) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_GE(st.models.size(), 1u);
    gemmi::Model& model = st.models[0];

    // Find two atoms in the same residue (try residue 10 first, then 1)
    const gemmi::Residue* res_ptr_10 = nullptr;
    const gemmi::Atom* atom1 = find_atom_with_residue(model, "A", 10, "N", &res_ptr_10);
    const gemmi::Residue* res_ptr_10_ca = nullptr;
    const gemmi::Atom* atom2 = find_atom_with_residue(model, "A", 10, "CA", &res_ptr_10_ca);

    if (!atom1 || !atom2) {
        res_ptr_10 = nullptr;
        res_ptr_10_ca = nullptr;
        atom1 = find_atom_with_residue(model, "A", 1, "N", &res_ptr_10);
        atom2 = find_atom_with_residue(model, "A", 1, "CA", &res_ptr_10_ca);
    }

    EXPECT_NE(atom1, nullptr);
    EXPECT_NE(atom2, nullptr);

    // In MMDB test, atom1->Het and atom2->Het were set to 1.
    // In gemmi, we pass atom1_is_het and atom2_is_het as parameters.
    bool atom1_is_het = true;
    bool atom2_is_het = true;

    // Setup protein geometry with dictionary
    coot::protein_geometry geom;
    geom.init_standard();
    bool have_dictionary = true;

    // Determine which residue seqnum we actually found
    int actual_seqnum = 10;
    if (res_ptr_10) {
        actual_seqnum = res_ptr_10->seqid.num.value;
    }

    // Case 1: two HET atoms in same residue (standard protein residue, likely no dictionary match)
    {
        gemmi::Residue* residue_1 = find_residue(model, "A", actual_seqnum);
        gemmi::Residue* residue_2 = residue_1; // same residue

        std::vector<std::pair<bool, gemmi::Residue *>> het_residues;
        int size_before = het_residues.size();
        bool result = add_bond_by_dictionary_maybe_gemmi(
            0, atom1, atom2, residue_1, residue_2, &het_residues, &geom, have_dictionary,
            atom1_is_het, atom2_is_het);
        int size_after = het_residues.size();
        EXPECT_EQ(result, true);
        EXPECT_EQ(size_before, 0);
        EXPECT_EQ(size_after, 1);
    }

    // Case 2: atoms in different residues (should short-circuit, return false)
    {
        const gemmi::Residue* res_ptr_other = nullptr;
        const gemmi::Atom* atom3 = find_atom_with_residue(model, "A", 11, "N", &res_ptr_other);
        if (!atom3) {
            res_ptr_other = nullptr;
            atom3 = find_atom_with_residue(model, "A", 2, "N", &res_ptr_other);
        }

        if (atom3) {
            int other_seqnum = res_ptr_other ? res_ptr_other->seqid.num.value : 2;
            gemmi::Residue* residue_1 = find_residue(model, "A", actual_seqnum);
            gemmi::Residue* residue_3 = find_residue(model, "A", other_seqnum);

            std::vector<std::pair<bool, gemmi::Residue *>> het_residues2;
            bool result = add_bond_by_dictionary_maybe_gemmi(
                0, atom1, atom3, residue_1, residue_3, &het_residues2, &geom, have_dictionary,
                atom1_is_het, atom2_is_het);
            EXPECT_EQ(result, false);
            EXPECT_EQ(het_residues2.size(), 0u);
        }
    }

    // Case 3: have_dictionary = false (should return false immediately)
    {
        gemmi::Residue* residue_1 = find_residue(model, "A", actual_seqnum);
        gemmi::Residue* residue_2 = residue_1; // same residue

        std::vector<std::pair<bool, gemmi::Residue *>> het_residues3;
        bool result = add_bond_by_dictionary_maybe_gemmi(
            0, atom1, atom2, residue_1, residue_2, &het_residues3, &geom, false,
            atom1_is_het, atom2_is_het);
        EXPECT_EQ(result, false);
        EXPECT_EQ(het_residues3.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}