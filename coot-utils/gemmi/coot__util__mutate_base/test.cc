#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__ideal_rna__get_standard_residue_instance/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__deep_copy_this_residue/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-dev/coot/utils/coot-utils.hh"
#include <iostream>

TEST(OracleTest, mutate_base) {
    // Load standard residues
    std::string standard_file = coot::package_data_dir() + "/standard-residues.pdb";
    gemmi::Structure st = gemmi::read_pdb_file(standard_file);
    ASSERT_FALSE(st.models.empty());
    const gemmi::Model &std_res_model = st.models[0];

    // Get standard adenine residue
    gemmi::Residue *std_A = coot::ideal_rna::get_standard_residue_instance_gemmi("A", std_res_model);
    ASSERT_NE(std_A, nullptr);

    // Get standard guanine residue
    gemmi::Residue *std_G = coot::ideal_rna::get_standard_residue_instance_gemmi("G", std_res_model);
    ASSERT_NE(std_G, nullptr);

    // Case 1: Mutate A to G (purine to purine)
    {
        gemmi::Residue *target_res = coot::util::deep_copy_this_residue_gemmi(std_A);
        gemmi::Residue *std_base = coot::util::deep_copy_this_residue_gemmi(std_G);

        int n_atoms_before = static_cast<int>(target_res->atoms.size());

        EXPECT_EQ(n_atoms_before, 22);
        EXPECT_EQ(target_res->name, "A");
        EXPECT_EQ(std_base->name, "G");

        coot::util::mutate_base_gemmi(*target_res, *std_base, true, false, 1.0f);

        int n_atoms_after = static_cast<int>(target_res->atoms.size());
        EXPECT_EQ(n_atoms_after, 23);
        EXPECT_EQ(target_res->name, "G");

        delete target_res;
        delete std_base;
    }

    // Case 2: Mutate A to A (same base - trivial)
    {
        gemmi::Residue *target_res = coot::util::deep_copy_this_residue_gemmi(std_A);
        gemmi::Residue *std_base = coot::util::deep_copy_this_residue_gemmi(std_A);

        int n_atoms_before = static_cast<int>(target_res->atoms.size());

        EXPECT_EQ(n_atoms_before, 22);
        EXPECT_EQ(target_res->name, "A");
        EXPECT_EQ(std_base->name, "A");

        coot::util::mutate_base_gemmi(*target_res, *std_base, true, false, 1.0f);

        int n_atoms_after = static_cast<int>(target_res->atoms.size());
        EXPECT_EQ(n_atoms_after, 22);
        EXPECT_EQ(target_res->name, "A");

        delete target_res;
        delete std_base;
    }

    // Case 3: Mutate G to C (purine to pyrimidine)
    {
        gemmi::Residue *target_res = coot::util::deep_copy_this_residue_gemmi(std_G);
        gemmi::Residue *std_base = coot::ideal_rna::get_standard_residue_instance_gemmi("C", std_res_model);
        ASSERT_NE(std_base, nullptr);

        int n_atoms_before = static_cast<int>(target_res->atoms.size());

        EXPECT_EQ(n_atoms_before, 23);
        EXPECT_EQ(target_res->name, "G");
        EXPECT_EQ(std_base->name, "C");

        coot::util::mutate_base_gemmi(*target_res, *std_base, true, false, 1.0f);

        int n_atoms_after = static_cast<int>(target_res->atoms.size());
        EXPECT_EQ(n_atoms_after, 20);
        EXPECT_EQ(target_res->name, "C");

        delete target_res;
        // std_base is NOT a deep_copy here (direct from get_standard_residue_instance), so don't delete
    }

    // Clean up the originals
    delete std_A;
    delete std_G;
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}