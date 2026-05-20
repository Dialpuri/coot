#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, dreiding_torsion_energy) {
    // Load PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());

    // Get ALA residue 21, chain A
    gemmi::Model& model = st.models[0];
    bool found = false;
    const gemmi::Atom *atom_N = nullptr;
    const gemmi::Atom *atom_CA = nullptr;
    const gemmi::Atom *atom_CB = nullptr;
    const gemmi::Atom *atom_C = nullptr;

    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == "A") {
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == 21) {
                    for (gemmi::Atom& atom : res.atoms) {
                        if (atom.name == "N") atom_N = &atom;
                        if (atom.name == "CA") atom_CA = &atom;
                        if (atom.name == "CB") atom_CB = &atom;
                        if (atom.name == "C") atom_C = &atom;
                    }
                    found = true;
                    break;
                }
            }
            if (found) break;
        }
    }

    EXPECT_NE(atom_N, nullptr);
    EXPECT_NE(atom_CA, nullptr);
    EXPECT_NE(atom_CB, nullptr);
    EXPECT_NE(atom_C, nullptr);

    // case: valid torsion (N-CA-CB-C in ALA residue 21, chain A)
    double energy = coot::protein_geometry::dreiding_torsion_energy_gemmi("ALA", 0, atom_N, atom_CA, atom_CB, atom_C);
    EXPECT_EQ(energy, 0);

    // case: residue code not in dictionary → returns 0 (indx == -1)
    double energy_unknown = coot::protein_geometry::dreiding_torsion_energy_gemmi("XYZUNKNOWN", 0, atom_N, atom_CA, atom_CB, atom_C);
    EXPECT_EQ(energy_unknown, 0);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}