#include <gtest/gtest.h>
#include <gemmi/cif.hpp>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, comp_tree_basic) {
    gemmi::cif::Loop loop;
    loop.tags.push_back("comp_id");
    loop.tags.push_back("atom_id");
    loop.tags.push_back("atom_back");
    loop.tags.push_back("atom_forward");
    loop.tags.push_back("connect_type");
    
    loop.values.push_back("LZA");
    loop.values.push_back("C1");
    loop.values.push_back("");
    loop.values.push_back("");
    loop.values.push_back("");
    
    loop.values.push_back("LZA");
    loop.values.push_back("C2");
    loop.values.push_back("C1");
    loop.values.push_back("");
    loop.values.push_back("covalent");
    
    int imol_enc = 0;
    EXPECT_NO_THROW(coot::protein_geometry::comp_tree_gemmi(loop, imol_enc));
}

TEST(OracleTest, comp_tree_empty_loop) {
    gemmi::cif::Loop loop;
    int imol_enc = 0;
    EXPECT_NO_THROW(coot::protein_geometry::comp_tree_gemmi(loop, imol_enc));
}

TEST(OracleTest, comp_tree_missing_columns) {
    gemmi::cif::Loop loop;
    loop.tags.push_back("comp_id");
    loop.tags.push_back("atom_id");
    
    loop.values.push_back("LZA");
    loop.values.push_back("C1");
    
    int imol_enc = 0;
    EXPECT_NO_THROW(coot::protein_geometry::comp_tree_gemmi(loop, imol_enc));
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}