#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: find an atom by chain_id, residue_seqid, atom_name (unpadded)
// Also returns a pointer to the parent residue.
static gemmi::Atom* find_atom_with_residue(
    gemmi::Structure &st,
    char chain_id,
    int res_seq,
    const std::string &atom_name,
    gemmi::Residue *&out_res)
{
    out_res = nullptr;
    if (st.models.empty()) return nullptr;
    gemmi::Model &mod = st.models[0];
    for (gemmi::Chain &ch : mod.chains) {
        if (ch.name != std::string(1, chain_id)) continue;
        for (gemmi::Residue &res : ch.residues) {
            if (res.seqid.num.value != res_seq) continue;
            for (gemmi::Atom &at : res.atoms) {
                if (at.name == atom_name) {
                    out_res = &res;
                    return &at;
                }
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, is_main_chain_p) {
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    // Case 1: CA atom (main chain) via atom pointer overload
    {
        gemmi::Residue *res = nullptr;
        gemmi::Atom *atom = find_atom_with_residue(st, 'A', 1, "CA", res);
        EXPECT_NE(atom, nullptr);
        bool result = coot::is_main_chain_p_gemmi(*atom, *res);
        EXPECT_TRUE(result);
    }

    // Case 2: CB atom (sidechain) via atom pointer overload
    {
        gemmi::Residue *res = nullptr;
        gemmi::Atom *atom = find_atom_with_residue(st, 'A', 1, "CB", res);
        EXPECT_NE(atom, nullptr);
        bool result = coot::is_main_chain_p_gemmi(*atom, *res);
        EXPECT_FALSE(result);
    }

    // Case 3: N atom (main chain) via atom pointer overload
    {
        gemmi::Residue *res = nullptr;
        gemmi::Atom *atom = find_atom_with_residue(st, 'A', 1, "N", res);
        EXPECT_NE(atom, nullptr);
        bool result = coot::is_main_chain_p_gemmi(*atom, *res);
        EXPECT_TRUE(result);
    }

    // Case 4: GLY HA2 (main chain for GLY) via atom pointer overload
    {
        gemmi::Residue *res = nullptr;
        gemmi::Atom *atom = find_atom_with_residue(st, 'A', 11, "HA2", res);
        EXPECT_NE(atom, nullptr);
        bool result = coot::is_main_chain_p_gemmi(*atom, *res);
        EXPECT_TRUE(result);
    }

    // Case 5: string overload with " CA "
    {
        std::string name = " CA ";
        bool result = coot::is_main_chain_p_gemmi(name);
        EXPECT_TRUE(result);
    }

    // Case 6: string overload with " CB " (sidechain in atom overload, but main-chain in string overload)
    {
        std::string name = " CB ";
        bool result = coot::is_main_chain_p_gemmi(name);
        EXPECT_TRUE(result);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
