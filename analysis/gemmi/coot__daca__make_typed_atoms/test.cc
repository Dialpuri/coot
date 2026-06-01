#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, make_typed_atoms_gemmi) {

    // Load PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model* model_p = nullptr;
    if (!st.models.empty()) {
        model_p = &st.models[0];
    }

    // Create protein_geometry and init
    coot::protein_geometry geom;
    geom.init_standard();

    // Case 1: valid model
    {
        std::vector<coot::daca::typed_atom_gemmi> typed_atoms =
            coot::daca::make_typed_atoms_gemmi(model_p, geom);

        EXPECT_NE(model_p, nullptr);
        if (model_p) {
            EXPECT_EQ(model_p->chains.size(), 2);
        }
        EXPECT_EQ(typed_atoms.size(), 2107u);

        // First few entries
        EXPECT_GE(typed_atoms.size(), 5u);
        {
            gemmi::Atom* at = typed_atoms[0].atom;
            EXPECT_NE(at, nullptr);
            EXPECT_EQ(at->name, "N");
            EXPECT_EQ(typed_atoms[0].residue->name, "MET");
            EXPECT_EQ(typed_atoms[0].residue->seqid.num.value, 1);
            EXPECT_EQ(typed_atoms[0].type, "NH1");
        }
        {
            gemmi::Atom* at = typed_atoms[1].atom;
            EXPECT_NE(at, nullptr);
            EXPECT_EQ(at->name, "CA");
            EXPECT_EQ(typed_atoms[1].residue->name, "MET");
            EXPECT_EQ(typed_atoms[1].residue->seqid.num.value, 1);
            EXPECT_EQ(typed_atoms[1].type, "CH1");
        }
        {
            gemmi::Atom* at = typed_atoms[2].atom;
            EXPECT_NE(at, nullptr);
            EXPECT_EQ(at->name, "C");
            EXPECT_EQ(typed_atoms[2].residue->name, "MET");
            EXPECT_EQ(typed_atoms[2].residue->seqid.num.value, 1);
            EXPECT_EQ(typed_atoms[2].type, "C");
        }
        {
            gemmi::Atom* at = typed_atoms[3].atom;
            EXPECT_NE(at, nullptr);
            EXPECT_EQ(at->name, "O");
            EXPECT_EQ(typed_atoms[3].residue->name, "MET");
            EXPECT_EQ(typed_atoms[3].residue->seqid.num.value, 1);
            EXPECT_EQ(typed_atoms[3].type, "O");
        }
        {
            gemmi::Atom* at = typed_atoms[4].atom;
            EXPECT_NE(at, nullptr);
            EXPECT_EQ(at->name, "CB");
            EXPECT_EQ(typed_atoms[4].residue->name, "MET");
            EXPECT_EQ(typed_atoms[4].residue->seqid.num.value, 1);
            EXPECT_EQ(typed_atoms[4].type, "CH2");
        }
    }

    // Case 2: null model — should return empty vector
    {
        std::vector<coot::daca::typed_atom_gemmi> typed_atoms =
            coot::daca::make_typed_atoms_gemmi(nullptr, geom);

        EXPECT_EQ(typed_atoms.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}