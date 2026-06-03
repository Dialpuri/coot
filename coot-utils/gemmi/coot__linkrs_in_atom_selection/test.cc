#include <gtest/gtest.h>
#include <gemmi/mmread.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, linkrs_in_atom_selection) {
    // Setup: load structure
    gemmi::Structure st = gemmi::read_structure_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");
    EXPECT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];

    coot::protein_geometry geom;

    EXPECT_EQ(st.connections.size(), 0u);

    // Case 1: Select atoms from LZA residue 1299, chain A (ligand of interest)
    {
        std::vector<gemmi::CRA> selected_atoms;
        for (auto& chain : model.chains) {
            if (chain.name != "A") continue;
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == 1299) {
                    for (auto& atom : res.atoms) {
                        selected_atoms.push_back(gemmi::CRA{&chain, &res, &atom});
                    }
                }
            }
        }

        EXPECT_EQ(selected_atoms.size(), 25u);

        bonded_pair_container_gemmi bpc = coot::linkrs_in_atom_selection_gemmi(
            st, selected_atoms, &geom);

        EXPECT_EQ(bpc.size(), 0u);
    }

    // Case 2: Empty selection — should return empty container
    {
        std::vector<gemmi::CRA> empty_atoms;
        bonded_pair_container_gemmi bpc = coot::linkrs_in_atom_selection_gemmi(
            st, empty_atoms, &geom);

        EXPECT_EQ(bpc.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}