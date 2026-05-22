#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include "function.hh"

// Test for is_inside_another_ligand_atom_gemmi
TEST(OracleTest, is_inside_another_ligand_atom) {
    // Create atom_overlaps_container_t using default constructor
    coot::atom_overlaps_container_t aoc;
    
    // Create test structure with gemmi
    gemmi::Structure st;
    gemmi::Model model;
    
    // Create chain and residue
    gemmi::Chain chain("A");
    gemmi::Residue res;
    res.seqid.num.value = 1;
    res.seqid.icode = ' ';
    res.name = "LIG";
    
    // Atom 1 at origin with radius 2.0
    gemmi::Atom atom1;
    atom1.name = "C1";
    atom1.pos = gemmi::Position(0.0, 0.0, 0.0);
    atom1.element = gemmi::Element("C");
    
    // Atom 2 at (3, 0, 0) with radius 1.5
    gemmi::Atom atom2;
    atom2.name = "C2";
    atom2.pos = gemmi::Position(3.0, 0.0, 0.0);
    atom2.element = gemmi::Element("C");
    
    res.atoms.push_back(atom1);
    res.atoms.push_back(atom2);
    chain.residues.push_back(res);
    model.chains.push_back(chain);
    st.models.push_back(model);
    
    aoc.model = &st.models[0];
    
    // Populate ligand_atom_neighbour_map
    // idx 0: neighbors of atom1
    std::vector<std::pair<gemmi::Atom *, double>> neighbors_0;
    neighbors_0.push_back(std::make_pair(&st.models[0].chains[0].residues[0].atoms[0], 2.0));
    neighbors_0.push_back(std::make_pair(&st.models[0].chains[0].residues[0].atoms[1], 1.5));
    aoc.ligand_atom_neighbour_map[0] = neighbors_0;
    
    // idx 1: neighbors of atom2
    std::vector<std::pair<gemmi::Atom *, double>> neighbors_1;
    neighbors_1.push_back(std::make_pair(&st.models[0].chains[0].residues[0].atoms[1], 1.5));
    aoc.ligand_atom_neighbour_map[1] = neighbors_1;
    
    // Test case 1: Point at origin (0,0,0) - inside atom1 (radius 2.0)
    {
        clipper::Coord_orth pt(0.0, 0.0, 0.0);
        bool result = aoc.is_inside_another_ligand_atom_gemmi(0, pt);
        EXPECT_EQ(result, true);
    }
    
    // Test case 2: Point at (5,0,0) - outside both atoms
    {
        clipper::Coord_orth pt(5.0, 0.0, 0.0);
        bool result = aoc.is_inside_another_ligand_atom_gemmi(0, pt);
        EXPECT_EQ(result, false);
    }
    
    // Test case 3: Point at (3.5,0,0) - inside atom2 (radius 1.5)
    {
        clipper::Coord_orth pt(3.5, 0.0, 0.0);
        bool result = aoc.is_inside_another_ligand_atom_gemmi(0, pt);
        EXPECT_EQ(result, true);
    }
    
    // Test case 4: Invalid index (negative)
    {
        clipper::Coord_orth pt(0.0, 0.0, 0.0);
        bool result = aoc.is_inside_another_ligand_atom_gemmi(-1, pt);
        EXPECT_EQ(result, false);
    }
    
    // Test case 5: Index not in map (999)
    {
        clipper::Coord_orth pt(0.0, 0.0, 0.0);
        bool result = aoc.is_inside_another_ligand_atom_gemmi(999, pt);
        EXPECT_EQ(result, false);
    }
    
    // Test case 6: Point at (5,0,0) for idx 1 - outside atom2
    {
        clipper::Coord_orth pt(5.0, 0.0, 0.0);
        bool result = aoc.is_inside_another_ligand_atom_gemmi(1, pt);
        EXPECT_EQ(result, false);
    }
    
    // Test case 7: Point at (3.2,0,0) for idx 1 - inside atom2 (radius 1.5)
    {
        clipper::Coord_orth pt(3.2, 0.0, 0.0);
        bool result = aoc.is_inside_another_ligand_atom_gemmi(1, pt);
        EXPECT_EQ(result, true);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}