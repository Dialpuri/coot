#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include <iostream>
#include <sstream>
#include <vector>

TEST(OracleTest, minimol_molecule_constructor) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];

    // Select only atoms from residue 1 of chain A
    std::vector<gemmi::CRA> cra_vec;
    std::vector<gemmi::Vec3> positions;

    for (auto& chain : model.chains) {
        if (chain.name != "A") continue;
        for (auto& res : chain.residues) {
            if (res.seqid.num.value != 1) continue;
            for (auto& atom : res.atoms) {
                cra_vec.push_back(gemmi::CRA{&chain, &res, &atom});
                positions.push_back(atom.pos);
            }
        }
    }

    int n_sel = static_cast<int>(cra_vec.size());
    EXPECT_EQ(n_sel, 8);
    EXPECT_GT(n_sel, 0);
    EXPECT_EQ(static_cast<int>(positions.size()), 8);

    // Case 1: Valid constructor call with matching sizes
    {
        coot::minimol::molecule_gemmi min_mol(cra_vec, n_sel, positions);

        EXPECT_EQ(static_cast<int>(min_mol.fragments.size()), 1);
        EXPECT_EQ(min_mol.have_cell, 0);
        EXPECT_EQ(min_mol.have_spacegroup, 0);

        EXPECT_EQ(min_mol.fragments[0].fragment_id, "A");
        EXPECT_EQ(static_cast<int>(min_mol.fragments[0].residues.size()), 2);

        // residue_0
        EXPECT_EQ(min_mol.fragments[0].residues[0].name, "");
        EXPECT_EQ(min_mol.fragments[0].residues[0].seqnum, 1);
        EXPECT_EQ(static_cast<int>(min_mol.fragments[0].residues[0].atoms.size()), 0);

        // residue_1
        EXPECT_EQ(min_mol.fragments[0].residues[1].name, "MET");
        EXPECT_EQ(min_mol.fragments[0].residues[1].seqnum, 1);
        EXPECT_EQ(static_cast<int>(min_mol.fragments[0].residues[1].atoms.size()), 8);
    }

    // Case 2: Mismatched sizes (atoms.size() != n_atoms) — triggers error guard
    {
        coot::minimol::molecule_gemmi min_mol2(cra_vec, 999, positions);
        EXPECT_EQ(static_cast<int>(min_mol2.fragments.size()), 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
