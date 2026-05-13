#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper to get atom_selection_container_t (defined in function.hh)

TEST(OracleTest, CellTranslationTests) {
    // Load test PDB
    std::string pdb_path = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";
    coot::atom_selection_container_t asc = coot::get_atom_selection(pdb_path, false, true, false);
    
    ASSERT_GT(asc.n_selected_atoms, 0) << "Failed to load PDB or no atoms selected";
    
    // Test 1: Point at (0, 0, 0) - origin should give (0,0,0)
    {
        coot::Cartesian point(0.0, 0.0, 0.0);
        coot::molecule_extents_t extents(asc, 0.0);
        coot::Cell_Translation trans = extents.coord_to_unit_cell_translations_gemmi(point, asc);
        EXPECT_EQ(trans.us, 0);
        EXPECT_EQ(trans.vs, 0);
        EXPECT_EQ(trans.ws, 0);
    }
    
    // Test 2: Point at actual atom coordinates - should map to appropriate unit cell
    {
        const auto& first_atom = asc.model->chains[0].residues[0].atoms[0];
        float x = first_atom.pos.x;
        float y = first_atom.pos.y;
        float z = first_atom.pos.z;
        coot::Cartesian point(x, y, z);
        coot::molecule_extents_t extents(asc, 0.0);
        coot::Cell_Translation trans = extents.coord_to_unit_cell_translations_gemmi(point, asc);
        // Based on oracle output: (0, 0, 1)
        EXPECT_EQ(trans.us, 0);
        EXPECT_EQ(trans.vs, 0);
        EXPECT_EQ(trans.ws, 1);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}