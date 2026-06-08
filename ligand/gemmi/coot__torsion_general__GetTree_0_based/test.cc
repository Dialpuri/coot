#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include <vector>
#include <string>
#include <cstdlib>

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    int r = RUN_ALL_TESTS();
    std::fflush(nullptr);
    _exit(r);
}

TEST(OracleTest, torsion_general_GetTree_0_based) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_EQ(st.models.size(), 1u);

    gemmi::Model& model = st.models[0];

    // Case 1: valid residue //A/1 with 4 torsion atoms N-CA-CB-CG
    {
        // Find residue A/1
        gemmi::Residue* res = nullptr;
        for (auto& chain : model.chains) {
            if (chain.name == "A") {
                for (auto& residue : chain.residues) {
                    if (residue.seqid.num.value == 1 && residue.name == "MET") {
                        res = &residue;
                        break;
                    }
                }
                break;
            }
        }
        EXPECT_NE(res, nullptr);

        // Find the 4 atoms N, CA, CB, CG
        gemmi::Atom* n_atom  = nullptr;
        gemmi::Atom* ca_atom = nullptr;
        gemmi::Atom* cb_atom = nullptr;
        gemmi::Atom* cg_atom = nullptr;

        for (auto& atom : res->atoms) {
            if (atom.name == "N"  ) n_atom  = &atom;
            if (atom.name == "CA" ) ca_atom = &atom;
            if (atom.name == "CB" ) cb_atom = &atom;
            if (atom.name == "CG" ) cg_atom = &atom;
        }

        // Oracle found all 4 atoms and built tree successfully
        EXPECT_NE(n_atom, nullptr);
        EXPECT_NE(ca_atom, nullptr);
        EXPECT_NE(cb_atom, nullptr);
        EXPECT_NE(cg_atom, nullptr);

        Tree tree = coot::torsion_general::GetTree_0_based_gemmi(res);
        EXPECT_EQ(tree.GetNumberOfVertices(), 8);
    }

    // Case 2: GLY //A/11 — only N,CA available (CB,CG missing)
    {
        // Find residue A/11
        gemmi::Residue* res = nullptr;
        for (auto& chain : model.chains) {
            if (chain.name == "A") {
                for (auto& residue : chain.residues) {
                    if (residue.seqid.num.value == 11 && residue.name == "GLY") {
                        res = &residue;
                        break;
                    }
                }
                break;
            }
        }
        EXPECT_NE(res, nullptr);

        // Find the atoms
        gemmi::Atom* n_atom  = nullptr;
        gemmi::Atom* ca_atom = nullptr;
        gemmi::Atom* cb_atom = nullptr;
        gemmi::Atom* cg_atom = nullptr;

        for (auto& atom : res->atoms) {
            if (atom.name == "N"  ) n_atom  = &atom;
            if (atom.name == "CA" ) ca_atom = &atom;
            if (atom.name == "CB" ) cb_atom = &atom;
            if (atom.name == "CG" ) cg_atom = &atom;
        }

        // Oracle: N and CA found, CB and CG not found
        EXPECT_NE(n_atom, nullptr);
        EXPECT_NE(ca_atom, nullptr);
        EXPECT_EQ(cb_atom, nullptr);
        EXPECT_EQ(cg_atom, nullptr);

        // With nullptr residue, tree is empty
        Tree tree = coot::torsion_general::GetTree_0_based_gemmi(nullptr);
        EXPECT_EQ(tree.GetNumberOfVertices(), 0);
    }
}
