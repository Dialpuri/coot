#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, mol_from_dictionary) {
    coot::protein_geometry geom;
    geom.init_standard();

    int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;

    // Case 1: valid standard amino acid (ALA)
    {
        std::string comp_id = "ALA";
        bool idealised = true;

        auto mol = coot::mol_from_dictionary_gemmi(geom, comp_id, imol_enc, idealised);
        EXPECT_NE(mol, nullptr);
        if (mol) {
            int n_models = mol->models.size();
            EXPECT_EQ(n_models, 1);
            if (n_models > 0) {
                gemmi::Model *model = &mol->models[0];
                if (model) {
                    int n_chains = model->chains.size();
                    EXPECT_EQ(n_chains, 1);
                    if (n_chains > 0) {
                        gemmi::Chain *chain = &model->chains[0];
                        if (chain) {
                            int n_res = chain->residues.size();
                            EXPECT_EQ(n_res, 1);
                            if (n_res > 0) {
                                gemmi::Residue *res = &chain->residues[0];
                                if (res) {
                                    int n_atoms = res->atoms.size();
                                    EXPECT_EQ(n_atoms, 13);
                                    EXPECT_STREQ(res->name.c_str(), "ALA");
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Case 2: invalid residue code — oracle returned non-null
    {
        std::string comp_id = "ZZZ";
        bool idealised = true;

        auto mol = coot::mol_from_dictionary_gemmi(geom, comp_id, imol_enc, idealised);
        EXPECT_NE(mol, nullptr);
        if (mol) {
            // no delete needed — unique_ptr handles cleanup
        }
    }

    // Case 3: non-idealised standard residue
    {
        std::string comp_id = "GLY";
        bool idealised = false;

        auto mol = coot::mol_from_dictionary_gemmi(geom, comp_id, imol_enc, idealised);
        EXPECT_NE(mol, nullptr);
        if (mol) {
            int n_models = mol->models.size();
            EXPECT_EQ(n_models, 1);
            if (n_models > 0) {
                gemmi::Model *model = &mol->models[0];
                if (model) {
                    int n_chains = model->chains.size();
                    EXPECT_EQ(n_chains, 1);
                    if (n_chains > 0) {
                        gemmi::Chain *chain = &model->chains[0];
                        if (chain) {
                            int n_res = chain->residues.size();
                            EXPECT_EQ(n_res, 1);
                            if (n_res > 0) {
                                gemmi::Residue *res = &chain->residues[0];
                                if (res) {
                                    int n_atoms = res->atoms.size();
                                    EXPECT_EQ(n_atoms, 10);
                                    EXPECT_STREQ(res->name.c_str(), "GLY");
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}