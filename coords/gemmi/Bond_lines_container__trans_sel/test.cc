#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, Bond_lines_container_trans_sel_gemmi) {
    // Setup: load PDB and collect atoms
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());

    // Collect all atoms (equivalent of make_asc)
    std::vector<gemmi::Atom> all_atoms;
    for (gemmi::Model& model : st.models) {
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& res : chain.residues) {
                for (gemmi::Atom& atom : res.atoms) {
                    all_atoms.push_back(atom);
                }
            }
        }
    }

    // Original: n_selected_atoms == 2107
    EXPECT_EQ(all_atoms.size(), 2107);

    // Case 1: identity symmetry transformation (isym=0, x=0, y=0, z=0; us=0, vs=0, ws=0)
    {
        std::vector<gemmi::Atom> trans_selection =
            coot::trans_sel_gemmi(&st, all_atoms, 0, 0, 0, 0, 0, 0, 0);

        ASSERT_FALSE(trans_selection.empty());
        EXPECT_EQ(trans_selection.size(), all_atoms.size());

        // Original: first->GetAtomName() == " N  "
        // gemmi Atom::name is unpadded ("N") vs MMDB's padded (" N  ")
        EXPECT_EQ(trans_selection[0].name, "N");
        // Original expected coords:
        EXPECT_NEAR(trans_selection[0].pos.x, 12.334, 1e-3);
        EXPECT_NEAR(trans_selection[0].pos.y, 2.772, 1e-3);
        EXPECT_NEAR(trans_selection[0].pos.z, 74.951, 1e-3);
    }

    // Case 2: non-identity symmetry (isym=1, x=1, y=0, z=0; us=-1, vs=0, ws=0)
    {
        std::vector<gemmi::Atom> trans_selection =
            coot::trans_sel_gemmi(&st, all_atoms, 1, 1, 0, 0, -1, 0, 0);

        ASSERT_FALSE(trans_selection.empty());
        EXPECT_EQ(trans_selection.size(), all_atoms.size());

        // Original: first->GetAtomName() == " N  "
        EXPECT_EQ(trans_selection[0].name, "N");
        // Original expected coords:
        EXPECT_NEAR(trans_selection[0].pos.x, 14.422, 1e-3);
        EXPECT_NEAR(trans_selection[0].pos.y, -2.772, 1e-3);
        EXPECT_NEAR(trans_selection[0].pos.z, 111.001, 1e-3);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}