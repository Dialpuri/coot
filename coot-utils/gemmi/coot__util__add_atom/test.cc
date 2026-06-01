#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: find a residue by chain and residue sequence number in gemmi
static gemmi::Residue* find_residue(gemmi::Model& model, const char* chain_id, int seqnum) {
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == chain_id) {
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == seqnum) {
                    return &res;
                }
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, add_atom) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];

    // --- Case 1: Valid call — add an atom to residue //A/1 (MET) using N, CA, C ---
    {
        gemmi::Residue* res = find_residue(model, "A", 1);
        ASSERT_NE(res, nullptr);

        int atoms_before = static_cast<int>(res->atoms.size());

        bool result = coot::util::add_atom_gemmi(
            res,
            " N  ",  // atom_name_1
            " CA ",  // atom_name_2
            " C  ",  // atom_name_3
            "",                 // alt_conf
            1.5,                // length (Angstroms)
            109.5,              // angle (degrees)
            180.0,              // torsion (degrees)
            "X1",               // new atom name
            "C",                // new atom element
            1.0f,               // occupancy
            20.0f               // B-factor
        );

        int atoms_after = static_cast<int>(res->atoms.size());

        EXPECT_TRUE(result);
        EXPECT_EQ(atoms_before, 8);
        EXPECT_EQ(atoms_after, 9);
    }

    // --- Case 2: Invalid atoms — try with atoms that don't exist ---
    {
        gemmi::Residue* res = find_residue(model, "A", 1);
        ASSERT_NE(res, nullptr);

        int atoms_before = static_cast<int>(res->atoms.size());

        bool result = coot::util::add_atom_gemmi(
            res,
            "ZZZ ",  // doesn't exist
            " CA ",
            " C  ",
            "",
            1.5, 109.5, 180.0,
            "X2", "C", 1.0f, 20.0f
        );

        int atoms_after = static_cast<int>(res->atoms.size());

        EXPECT_FALSE(result);
        EXPECT_EQ(atoms_before, 9);
        EXPECT_EQ(atoms_after, 9);
    }

    // --- Case 3: Null residue ---
    {
        bool result = coot::util::add_atom_gemmi(
            nullptr,
            " N  ", " CA ", " C  ",
            "",
            1.5, 109.5, 180.0,
            "X3", "C", 1.0f, 20.0f
        );

        EXPECT_FALSE(result);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}