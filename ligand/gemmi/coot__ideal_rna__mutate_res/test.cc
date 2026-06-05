#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: create a residue with one atom and add to a model chain
static void add_res_to_model(gemmi::Model &model, const char *name) {
    gemmi::Chain chain;
    chain.name = "X";
    gemmi::Residue res;
    res.name = name;
    res.seqid.num.value = 1;
    gemmi::Atom atom;
    atom.name = "C1'";
    atom.pos = gemmi::Position(0, 0, 0);
    res.atoms.push_back(std::move(atom));
    chain.residues.push_back(std::move(res));
    model.chains.push_back(std::move(chain));
}

TEST(OracleTest, mutate_res) {
    // Create a standard_residues model with nucleic acid residues
    gemmi::Model standard_residues;

    // Add G, A, C, U, DA, DG, DC, DT residues
    add_res_to_model(standard_residues, "G");
    add_res_to_model(standard_residues, "A");
    add_res_to_model(standard_residues, "C");
    add_res_to_model(standard_residues, "U");
    add_res_to_model(standard_residues, "DA");
    add_res_to_model(standard_residues, "DG");
    add_res_to_model(standard_residues, "DC");
    add_res_to_model(standard_residues, "DT");

    // Verify we can find residues
    gemmi::Residue *std_g = coot::ideal_rna::get_standard_residue_instance_gemmi("G", standard_residues);
    ASSERT_NE(std_g, nullptr);
    delete std_g;

    // Test case 1: RNA 'g' -> G (should succeed, status=1)
    {
        gemmi::Residue *res = coot::ideal_rna::get_standard_residue_instance_gemmi("G", standard_residues);
        ASSERT_NE(res, nullptr);
        res->seqid.num.value = 1;
        res->name = "G";

        char base = 'g';
        bool is_dna_flag = false;

        int status = coot::ideal_rna::mutate_res_gemmi(res, base, is_dna_flag, standard_residues);

        EXPECT_EQ(status, 1);
        EXPECT_STREQ(res->name.c_str(), "G");

        delete res;
    }

    // Test case 2: DNA 't' -> DT (should succeed, status=1)
    {
        gemmi::Residue *res = coot::ideal_rna::get_standard_residue_instance_gemmi("G", standard_residues);
        ASSERT_NE(res, nullptr);
        res->seqid.num.value = 2;
        res->name = "G";

        char base = 't';
        bool is_dna_flag = true;

        int status = coot::ideal_rna::mutate_res_gemmi(res, base, is_dna_flag, standard_residues);

        EXPECT_EQ(status, 1);
        EXPECT_STREQ(res->name.c_str(), "G");

        delete res;
    }

    // Test case 3: invalid base 'x' -> should fail (status=0)
    {
        gemmi::Residue *res = coot::ideal_rna::get_standard_residue_instance_gemmi("G", standard_residues);
        ASSERT_NE(res, nullptr);
        res->seqid.num.value = 3;
        res->name = "G";

        char base = 'x';
        bool is_dna_flag = false;

        int status = coot::ideal_rna::mutate_res_gemmi(res, base, is_dna_flag, standard_residues);

        EXPECT_EQ(status, 0);
        EXPECT_STREQ(res->name.c_str(), "G");

        delete res;
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}