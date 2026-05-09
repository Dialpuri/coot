#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, GetSVGFor2DLigandEnvironmentView) {
    // Read the PDB file using gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Test 1: //A/10 (ILE, terminal residue, should fail due to OXT)
    {
        std::string cid = "//A/10";
        bool add_key = false;
        coot::molecule_t mol;
        std::string svg = mol.get_svg_for_2d_ligand_environment_view_gemmi(cid, nullptr, add_key, st);
        EXPECT_EQ(svg.length(), 0U);
    }

    // Test 2: //A/50 (ARG, internal, should fail due to OXT)
    {
        std::string cid = "//A/50";
        bool add_key = false;
        coot::molecule_t mol;
        std::string svg = mol.get_svg_for_2d_ligand_environment_view_gemmi(cid, nullptr, add_key, st);
        EXPECT_EQ(svg.length(), 0U);
    }

    // Test 3: //A/62 (ASN, internal, should fail due to OXT)
    {
        std::string cid = "//A/62";
        bool add_key = false;
        coot::molecule_t mol;
        std::string svg = mol.get_svg_for_2d_ligand_environment_view_gemmi(cid, nullptr, add_key, st);
        EXPECT_EQ(svg.length(), 0U);
    }

    // Test 4: //A/50 with add_key=true (should also fail)
    {
        std::string cid = "//A/50";
        bool add_key = true;
        coot::molecule_t mol;
        std::string svg = mol.get_svg_for_2d_ligand_environment_view_gemmi(cid, nullptr, add_key, st);
        EXPECT_EQ(svg.length(), 0U);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}