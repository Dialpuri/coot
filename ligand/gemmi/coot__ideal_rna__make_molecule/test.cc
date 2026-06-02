#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "coot/utils/coot-utils.hh"
#include "function.hh"

TEST(OracleTest, ideal_rna) {
    // Load standard residues
    gemmi::Structure st = gemmi::read_pdb_file(coot::package_data_dir() + "/standard-residues.pdb");
    ASSERT_EQ(st.models.size(), 1u);
    gemmi::Model& standard_residues = st.models[0];

    // Case 1: RNA, A-form, single-stranded, sequence "ACGU"
    {
        std::string seq = "ACGU";
        std::string RNA_or_DNA = "RNA";
        std::string form = "A";
        bool single_stranded = true;

        auto mol = coot::ideal_rna::make_molecule_gemmi(seq, RNA_or_DNA, form, single_stranded, standard_residues);

        EXPECT_NE(mol, nullptr);

        if (mol) {
            int n_models = static_cast<int>(mol->models.size());
            EXPECT_EQ(n_models, 1);
            if (n_models > 0) {
                gemmi::Model& mod = mol->models[0];
                int n_chains = static_cast<int>(mod.chains.size());
                EXPECT_EQ(n_chains, 0);
            }
        }
    }

    // Case 2: DNA, B-form, double-stranded, sequence "AA"
    {
        std::string seq = "AA";
        std::string RNA_or_DNA = "DNA";
        std::string form = "B";
        bool single_stranded = false;

        auto mol = coot::ideal_rna::make_molecule_gemmi(seq, RNA_or_DNA, form, single_stranded, standard_residues);

        EXPECT_NE(mol, nullptr);

        if (mol) {
            int n_models = static_cast<int>(mol->models.size());
            EXPECT_EQ(n_models, 1);
            if (n_models > 0) {
                gemmi::Model& mod = mol->models[0];
                int n_chains = static_cast<int>(mod.chains.size());
                EXPECT_EQ(n_chains, 0);
            }
        }
    }

    // Case 3: empty sequence — should return nullptr
    {
        std::string seq = "";
        std::string RNA_or_DNA = "RNA";
        std::string form = "A";
        bool single_stranded = true;

        auto mol = coot::ideal_rna::make_molecule_gemmi(seq, RNA_or_DNA, form, single_stranded, standard_residues);

        EXPECT_EQ(mol, nullptr);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}