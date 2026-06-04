#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Find residue in gemmi Structure by chain name and sequence number
static const gemmi::Residue* find_residue(const gemmi::Structure& st, const std::string& chain_name, int seq_num) {
    for (const auto& model : st.models) {
        for (const auto& chain : model.chains) {
            if (chain.name != chain_name) continue;
            for (const auto& res : chain.residues) {
                if (res.seqid.num.value == seq_num) return &res;
            }
        }
    }
    return nullptr;
}

// Build atom_torsion_data_t entries that mimic what link_by_torsion_t loads for "alpha1-4".
// The critical detail: prior_atom_1 references "C2" or "C1" or "O5" — sugar atoms that
// a protein residue (MET) does NOT have, so pos_gemmi will throw.
static std::vector<atom_torsion_data_t> make_alpha1_4_atoms() {
    // Minimal set of atoms for an NAG/MAN residue linked alpha1-4.
    // First atom references C2 on the base residue (sugar atom not in MET).
    std::vector<atom_torsion_data_t> atoms;
    atom_torsion_data_t a1;
    a1.atom_name = "C1";
    a1.element = "C";
    a1.prior_atom_1 = {true, "C2"};  // sugar atom — MET doesn't have this → throw
    a1.prior_atom_2 = {true, "C1"};
    a1.prior_atom_3 = {true, "O5"};
    a1.bond_length = 1.53;
    a1.angle_deg = 109.5;
    a1.torsion_deg = -60.0;
    atoms.push_back(a1);
    return atoms;
}

TEST(OracleTestGemmi, link_by_torsion) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Case 1: valid link type + residue type, but base residue is a protein (not sugar)
    //          -> throws runtime_error because base residue lacks sugar atoms (C2, C1, O5)
    {
        const gemmi::Residue* base_res = find_residue(st, "A", 1);
        EXPECT_NE(base_res, nullptr);
        if (base_res == nullptr) return;
        EXPECT_EQ(base_res->name, "MET");

        float b_factor = 20.0f;
        auto atoms = make_alpha1_4_atoms();
        EXPECT_THROW(coot::make_residue_gemmi("NAG", 999, atoms, *base_res, b_factor), std::runtime_error);
    }

    // Case 2: invalid link type (file not found) — geom_atom_torsions empty -> returns NULL
    {
        const gemmi::Residue* base_res = find_residue(st, "A", 1);
        EXPECT_NE(base_res, nullptr);
        if (base_res == nullptr) return;

        std::vector<atom_torsion_data_t> atoms; // empty — mimics !filled()
        gemmi::Residue* result = coot::make_residue_gemmi("FOO", 1, atoms, *base_res, 20.0f);
        EXPECT_EQ(result, nullptr);
    }

    // Case 3: valid link type with different sugar residue type (MAN)
    {
        const gemmi::Residue* base_res = find_residue(st, "A", 1);
        EXPECT_NE(base_res, nullptr);
        if (base_res == nullptr) return;
        EXPECT_EQ(base_res->name, "MET");

        float b_factor = 30.0f;
        auto atoms = make_alpha1_4_atoms();
        EXPECT_THROW(coot::make_residue_gemmi("MAN", 999, atoms, *base_res, b_factor), std::runtime_error);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}