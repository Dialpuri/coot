#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "/lmb/home/jdialpuri/Development/coot-dev/coot/utils/coot-utils.hh"

TEST(OracleTest, ideal_rna_get_standard_residue_instance) {
    // Load standard residues
    std::string standard_file = coot::package_data_dir() + "/standard-residues.pdb";

    gemmi::Structure st = gemmi::read_pdb_file(standard_file);
    ASSERT_FALSE(st.models.empty());
    const gemmi::Model& std_res_model = st.models[0];

    // Case 1: RNA residue "A" (adenine)
    {
        std::string residue_type_in = "A";
        gemmi::Residue *res = coot::ideal_rna::get_standard_residue_instance_gemmi(residue_type_in, std_res_model);
        EXPECT_NE(res, nullptr);
        if (res) {
            EXPECT_EQ(res->name, "A");
            EXPECT_EQ(res->atoms.size(), 22);
            delete res;
        }
    }

    // Case 2: DNA residue "DA" (deoxyadenine)
    {
        std::string residue_type_in = "DA";
        gemmi::Residue *res = coot::ideal_rna::get_standard_residue_instance_gemmi(residue_type_in, std_res_model);
        EXPECT_NE(res, nullptr);
        if (res) {
            EXPECT_EQ(res->name, "DA");
            EXPECT_EQ(res->atoms.size(), 21);
            delete res;
        }
    }

    // Case 3: Invalid residue type "INVALID"
    {
        std::string residue_type_in = "INVALID";
        gemmi::Residue *res = coot::ideal_rna::get_standard_residue_instance_gemmi(residue_type_in, std_res_model);
        EXPECT_EQ(res, nullptr);
        if (res) {
            delete res;
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}