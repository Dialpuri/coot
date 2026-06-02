#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// ... TEST(...) blocks here ...

TEST(OracleTest, mol_from_dictionary) {
    coot::protein_geometry geom;
    geom.init_standard();

    int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;

    // Case 1: Valid standard amino acid "ALA"
    {
        std::string comp_id = "ALA";
        bool idealised = true;
        gemmi::Structure *mol = coot::mol_from_dictionary_gemmi(geom, comp_id, imol_enc, idealised);
        EXPECT_NE(mol, nullptr);
        if (mol) {
            EXPECT_EQ(mol->models.size(), 1u);
            gemmi::Model &model = mol->models[0];
            if (true) {
                EXPECT_EQ(model.chains.size(), 1u);
            }
            delete mol;
        }
    }

    // Case 2: Another valid standard amino acid "GLY"
    {
        std::string comp_id = "GLY";
        bool idealised = true;
        gemmi::Structure *mol = coot::mol_from_dictionary_gemmi(geom, comp_id, imol_enc, idealised);
        EXPECT_NE(mol, nullptr);
        if (mol) {
            delete mol;
        }
    }

    // Case 3: Invalid/non-existent comp_id
    {
        std::string comp_id = "ZZZ";
        bool idealised = true;
        gemmi::Structure *mol = coot::mol_from_dictionary_gemmi(geom, comp_id, imol_enc, idealised);
        EXPECT_NE(mol, nullptr);
        if (mol) {
            delete mol;
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}