#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(AtomQuadTest, SetupChiralQuad) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    const gemmi::Residue* res1 = nullptr;
    const gemmi::Residue* res2 = nullptr;
    
    for (const auto& model : st.models) {
        for (const auto& chain : model.chains) {
            for (const auto& res : chain.residues) {
                if (res.seqid.num.value == 10 && chain.name == "A") {
                    res1 = &res;
                }
                if (res.seqid.num.value == 11 && chain.name == "A") {
                    res2 = &res;
                }
            }
        }
    }
    
    EXPECT_NE(res1, nullptr);
    EXPECT_NE(res2, nullptr);
    
    if (res1 && res2) {
        std::vector<std::string> chiral_atom_names;
        coot::atom_quad::atom_quad quad = coot::atom_quad::setup_chiral_quad_gemmi(res1, res2, "", chiral_atom_names);
        
        EXPECT_EQ(quad.atom_1, nullptr);
        EXPECT_EQ(quad.atom_2, nullptr);
        EXPECT_EQ(quad.atom_3, nullptr);
        EXPECT_EQ(quad.atom_4, nullptr);
        EXPECT_EQ(quad.name, "");
    }
    
    if (res1 && res2) {
        std::vector<std::string> chiral_atom_names;
        chiral_atom_names.push_back("C1");
        chiral_atom_names.push_back("C2");
        chiral_atom_names.push_back("C3");
        
        coot::atom_quad::atom_quad quad = coot::atom_quad::setup_chiral_quad_gemmi(res1, res2, "O1", chiral_atom_names);
        
        EXPECT_EQ(quad.atom_1, nullptr);
        EXPECT_EQ(quad.atom_2, nullptr);
        EXPECT_EQ(quad.atom_3, nullptr);
        EXPECT_EQ(quad.atom_4, nullptr);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}