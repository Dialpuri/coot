#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

using namespace coot;

TEST(OracleTest, dictionary_residue_restraints_t_make_graph_gemmi) {
    coot::dictionary_residue_restraints_t drr;

    // Set comp_id
    drr.residue_info.comp_id = "TEST";

    // Add atoms: C1, N1, C2
    {
        dict_atom atom;
        atom.atom_id_4c = "C1  ";
        atom.type_symbol = " C";
        drr.atom_info.push_back(atom);
    }
    {
        dict_atom atom;
        atom.atom_id_4c = "N1  ";
        atom.type_symbol = " N";
        drr.atom_info.push_back(atom);
    }
    {
        dict_atom atom;
        atom.atom_id_4c = "C2  ";
        atom.type_symbol = " C";
        drr.atom_info.push_back(atom);
    }

    // Add bonds: C1-N1, N1-C2
    {
        dict_bond_restraint_t brt("C1  ", "N1  ", "single");
        drr.bond_restraint.push_back(brt);
    }
    {
        dict_bond_restraint_t brt("N1  ", "C2  ", "single");
        drr.bond_restraint.push_back(brt);
    }

    // Case 1: make_graph without hydrogens
    {
        bool use_h = false;
        std::unique_ptr<dict_graph_t> g = coot::make_graph_gemmi(drr, use_h);
        EXPECT_NE(g, nullptr);
        if (g) {
            EXPECT_EQ(g->GetNofVertices(), 3);
            EXPECT_EQ(g->GetNofEdges(), 2);
        }
    }

    // Case 2: with hydrogens (same result since no H atoms)
    {
        bool use_h = true;
        std::unique_ptr<dict_graph_t> g = coot::make_graph_gemmi(drr, use_h);
        EXPECT_NE(g, nullptr);
        if (g) {
            EXPECT_EQ(g->GetNofVertices(), 3);
            EXPECT_EQ(g->GetNofEdges(), 2);
        }
    }

    // Case 3: add a hydrogen atom and test filtering
    {
        dict_atom atom;
        atom.atom_id_4c = "H1  ";
        atom.type_symbol = " H";
        drr.atom_info.push_back(atom);
    }
    {
        dict_bond_restraint_t brt("C2  ", "H1  ", "single");
        drr.bond_restraint.push_back(brt);
    }

    // Without hydrogens: should have 3 vertices, 2 edges
    {
        bool use_h = false;
        EXPECT_EQ(drr.atom_info.size(), 4);
        std::unique_ptr<dict_graph_t> g = coot::make_graph_gemmi(drr, use_h);
        EXPECT_NE(g, nullptr);
        if (g) {
            EXPECT_EQ(g->GetNofVertices(), 3);
            EXPECT_EQ(g->GetNofEdges(), 2);
        }
    }

    // With hydrogens: should have 4 vertices, 3 edges
    {
        bool use_h = true;
        std::unique_ptr<dict_graph_t> g = coot::make_graph_gemmi(drr, use_h);
        EXPECT_NE(g, nullptr);
        if (g) {
            EXPECT_EQ(g->GetNofVertices(), 4);
            EXPECT_EQ(g->GetNofEdges(), 3);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}