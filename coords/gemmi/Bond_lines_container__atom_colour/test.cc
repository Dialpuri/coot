#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, Bond_lines_container_atom_colour) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    gemmi::Model& model = st.models[0];

    // Find the CA atom of residue 10 in chain A (seqid 10, ILE)
    gemmi::CRA cra = {};
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name != "A") continue;
        for (gemmi::Residue& res : chain.residues) {
            if (res.seqid.num.value != 10) continue;
            for (gemmi::Atom& atom : res.atoms) {
                if (atom.name == "CA") {
                    cra = gemmi::CRA{&chain, &res, &atom};
                    goto found;
                }
            }
        }
    }
found:
    ASSERT_TRUE(cra.atom != nullptr);

    int udd_handle = -1;
    my_atom_colour_map_t *atom_colour_map_p = nullptr;

    // Case 1: COLOUR_BY_ATOM_TYPE for Carbon atom (CA)
    {
        int col = atom_colour_gemmi(cra, BC_ATOM_TYPE, udd_handle, atom_colour_map_p);
        EXPECT_EQ(col, 0);
    }

    // Case 2: COLOUR_BY_MOLECULE (returns 0 — one colour fits all)
    {
        int col = atom_colour_gemmi(cra, BC_MOLECULE, udd_handle, atom_colour_map_p);
        EXPECT_EQ(col, 0);
    }

    // Case 3: COLOUR_BY_B_FACTOR
    {
        int col = atom_colour_gemmi(cra, BC_B_FACTOR, udd_handle, atom_colour_map_p);
        EXPECT_EQ(col, 18);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}