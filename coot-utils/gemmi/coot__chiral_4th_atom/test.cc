#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, Chiral4thAtom) {
    // Load the test PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty()) << "Failed to load PDB file";

    // case 1: ILE 10 - chiral CB with CG1, CG2, CD1
    {
        bool found = false;
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                for (auto& res : chain.residues) {
                    if (res.name == "ILE" && res.seqid.num.value == 10) {
                        found = true;
                        
                        gemmi::Atom* at_cb = nullptr;
                        gemmi::Atom* at_cg1 = nullptr;
                        gemmi::Atom* at_cg2 = nullptr;
                        gemmi::Atom* at_cd1 = nullptr;
                        
                        for (auto& atom : res.atoms) {
                            if (atom.name == "CB") at_cb = &atom;
                            if (atom.name == "CG1") at_cg1 = &atom;
                            if (atom.name == "CG2") at_cg2 = &atom;
                            if (atom.name == "CD1") at_cd1 = &atom;
                        }
                        
                        ASSERT_NE(at_cb, nullptr);
                        ASSERT_NE(at_cg1, nullptr);
                        ASSERT_NE(at_cg2, nullptr);
                        ASSERT_NE(at_cd1, nullptr);
                        
                        gemmi::Atom* at_4th = coot::chiral_4th_atom_gemmi(&res, at_cb, at_cg1, at_cg2, at_cd1);
                        
                        EXPECT_EQ(at_4th, nullptr) << "Expected no fourth atom for ILE 10";
                    }
                }
            }
        }
        EXPECT_TRUE(found) << "ILE 10 not found";
    }

    // case 2: ARG 50 - chiral CB with CG, CD, NE
    {
        bool found = false;
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                for (auto& res : chain.residues) {
                    if (res.name == "ARG" && res.seqid.num.value == 50) {
                        found = true;
                        
                        gemmi::Atom* at_cb = nullptr;
                        gemmi::Atom* at_cg = nullptr;
                        gemmi::Atom* at_cd = nullptr;
                        gemmi::Atom* at_ne = nullptr;
                        
                        for (auto& atom : res.atoms) {
                            if (atom.name == "CB") at_cb = &atom;
                            if (atom.name == "CG") at_cg = &atom;
                            if (atom.name == "CD") at_cd = &atom;
                            if (atom.name == "NE") at_ne = &atom;
                        }
                        
                        ASSERT_NE(at_cb, nullptr);
                        ASSERT_NE(at_cg, nullptr);
                        ASSERT_NE(at_cd, nullptr);
                        ASSERT_NE(at_ne, nullptr);
                        
                        gemmi::Atom* at_4th = coot::chiral_4th_atom_gemmi(&res, at_cb, at_cg, at_cd, at_ne);
                        
                        EXPECT_EQ(at_4th, nullptr) << "Expected no fourth atom for ARG 50";
                    }
                }
            }
        }
        EXPECT_TRUE(found) << "ARG 50 not found";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}