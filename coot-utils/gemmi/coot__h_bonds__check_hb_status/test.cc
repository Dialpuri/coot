#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include <cstdlib>

TEST(OracleTest, check_hb_status) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());

    coot::protein_geometry geom;
    geom.init_standard();

    int imol = 0;
    gemmi::Model& model = st.models[0];

    // Case 1: select all atoms of residue with seqid 10 in the first chain
    // Original: mol->SelectAtoms(selhnd, 0, "*", 10, "*", 10, "*", "*", "*", "*", "*");
    {
        std::vector<gemmi::CRA> sel_atoms;
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == 10) {
                    for (gemmi::Atom& atom : res.atoms) {
                        sel_atoms.push_back({&chain, &res, &atom});
                    }
                }
            }
        }

        std::pair<bool, std::map<gemmi::Atom*, int>> status =
            coot::h_bonds::check_hb_status_gemmi(sel_atoms, geom, imol);
        EXPECT_TRUE(status.first);
        // In the gemmi version, status.second is the map (not a UDD handle).
        // The original asserted status.second == 16777218 (the UDD handle value).
        // Here we assert the map is non-empty, consistent with status.first == true.
        EXPECT_FALSE(status.second.empty());
    }

    // Case 2: select all atoms in the molecule
    // Original: mol->SelectAtoms(selhnd, 0, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES, "*", "*", "*", "*", "*");
    {
        std::vector<gemmi::CRA> sel_atoms;
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& res : chain.residues) {
                for (gemmi::Atom& atom : res.atoms) {
                    sel_atoms.push_back({&chain, &res, &atom});
                }
            }
        }

        std::pair<bool, std::map<gemmi::Atom*, int>> status =
            coot::h_bonds::check_hb_status_gemmi(sel_atoms, geom, imol);
        EXPECT_TRUE(status.first);
        // Original: EXPECT_EQ(status.second, 16777218) — UDD handle
        // Adapted: map is non-empty since atoms were found and marked
        EXPECT_FALSE(status.second.empty());
    }

    // Case 3: empty selection (no atoms) — edge case, should return false
    // Original: mol->SelectAtoms(selhnd, 0, "*", 99999, "*", 99999, "*", "*", "*", "*", "*");
    {
        std::vector<gemmi::CRA> sel_atoms;  // empty — residue 99999 does not exist

        std::pair<bool, std::map<gemmi::Atom*, int>> status =
            coot::h_bonds::check_hb_status_gemmi(sel_atoms, geom, imol);
        EXPECT_FALSE(status.first);
        // Original: EXPECT_EQ(status.second, 16777218) — UDD handle still returned even for empty selection
        // Adapted: map is empty because no atoms were selected
        EXPECT_TRUE(status.second.empty());
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}