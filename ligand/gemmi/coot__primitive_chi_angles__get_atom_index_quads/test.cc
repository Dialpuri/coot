#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

using coot::atom_name_quad;
using coot::atom_index_quad;
using coot::primitive_chi_angles::get_atom_index_quads_gemmi;

TEST(OracleTest, PrimitiveChiAnglesGetAtomIndexQuads) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    gemmi::Chain& chain = st.models[0].chains[0];
    gemmi::Residue& res = chain.residues[9];
    
    EXPECT_EQ(res.name, "ILE");
    
    std::vector<std::string> atom_names;
    for (const auto& atom : res.atoms) {
        atom_names.push_back(atom.name);
    }
    
    std::vector<atom_name_quad> atom_name_quads;
    atom_name_quads.push_back(atom_name_quad("N", "CA", "CB", "CG1"));
    atom_name_quads.push_back(atom_name_quad("CA", "CB", "CG1", "CD1"));
    
    EXPECT_EQ(atom_name_quads.size(), 2u);
    
    EXPECT_EQ(atom_name_quads[0].atom_name(0), "N");
    EXPECT_EQ(atom_name_quads[0].atom_name(1), "CA");
    EXPECT_EQ(atom_name_quads[0].atom_name(2), "CB");
    EXPECT_EQ(atom_name_quads[0].atom_name(3), "CG1");
    
    EXPECT_EQ(atom_name_quads[1].atom_name(0), "CA");
    EXPECT_EQ(atom_name_quads[1].atom_name(1), "CB");
    EXPECT_EQ(atom_name_quads[1].atom_name(2), "CG1");
    EXPECT_EQ(atom_name_quads[1].atom_name(3), "CD1");
    
    std::vector<atom_index_quad> result = get_atom_index_quads_gemmi(atom_name_quads, atom_names);
    
    EXPECT_EQ(result.size(), 2u);
    EXPECT_EQ(atom_name_quads.size(), 2u);
    
    if (result.size() > 0) {
        EXPECT_EQ(result[0].index1, 0);
        EXPECT_EQ(result[0].index2, 1);
        EXPECT_EQ(result[0].index3, 4);
        EXPECT_EQ(result[0].index4, 5);
    }
    
    if (result.size() > 1) {
        EXPECT_EQ(result[1].index1, 1);
        EXPECT_EQ(result[1].index2, 4);
        EXPECT_EQ(result[1].index3, 5);
        EXPECT_EQ(result[1].index4, 7);
    }
    
    {
        std::vector<atom_name_quad> empty_quads;
        std::vector<atom_index_quad> empty_result = get_atom_index_quads_gemmi(empty_quads, atom_names);
        EXPECT_EQ(empty_result.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}