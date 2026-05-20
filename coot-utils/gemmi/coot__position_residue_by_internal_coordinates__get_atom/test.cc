#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, position_residue_by_internal_coordinates_get_atom) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Get two residues - ASN at position 100 and ASN at position 101
    gemmi::Residue *res1 = nullptr;
    gemmi::Residue *res2 = nullptr;
    
    for (auto& model : st.models)
        for (auto& chain : model.chains)
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == 100 && chain.name == "A") {
                    res1 = &res;
                }
                if (res.seqid.num.value == 101 && chain.name == "A") {
                    res2 = &res;
                }
            }
    
    EXPECT_NE(res1, nullptr);
    if (res1) {
        EXPECT_EQ(res1->seqid.num.value, 100);
        EXPECT_EQ(res1->name, "PRO");
    }
    
    EXPECT_NE(res2, nullptr);
    if (res2) {
        EXPECT_EQ(res2->seqid.num.value, 101);
        EXPECT_EQ(res2->name, "LEU");
    }
    
    if (!res1 || !res2) {
        return;
    }
    
    // Create atom_name_quad with valid atom names
    coot::atom_name_quad quad("OD1", "ND2", "CB", "CG");
    quad.set_atom_residue_index(0, 1);
    quad.set_atom_residue_index(1, 1);
    quad.set_atom_residue_index(2, 2);
    quad.set_atom_residue_index(3, 2);
    
    EXPECT_EQ(quad.format(), "(OD1 ND2 CB CG)");
    EXPECT_EQ(quad.atom_residue_index[0], 1);
    EXPECT_EQ(quad.atom_residue_index[1], 1);
    EXPECT_EQ(quad.atom_residue_index[2], 2);
    EXPECT_EQ(quad.atom_residue_index[3], 2);
    
    // Test case 1: Call get_atom with atom_index=0
    {
        gemmi::Atom* at = coot::position_residue_by_internal_coordinates_get_atom_gemmi(res1, res2, quad, 0);
        EXPECT_EQ(at, nullptr);
    }
    
    // Test case 2: Call get_atom with atom_index=1
    {
        gemmi::Atom* at = coot::position_residue_by_internal_coordinates_get_atom_gemmi(res1, res2, quad, 1);
        EXPECT_EQ(at, nullptr);
    }
    
    // Test case 3: Call get_atom with atom_index=2 (from res2)
    {
        gemmi::Atom* at = coot::position_residue_by_internal_coordinates_get_atom_gemmi(res1, res2, quad, 2);
        EXPECT_EQ(at, nullptr);
    }
    
    // Test case 4: NULL residue (should return NULL)
    {
        gemmi::Residue* null_res = nullptr;
        gemmi::Atom* at = coot::position_residue_by_internal_coordinates_get_atom_gemmi(null_res, res2, quad, 0);
        EXPECT_EQ(at, nullptr);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}