#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, is_inside_another_ligand_atom) {
    // Load PDB to get atoms
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Get a model
    ASSERT_FALSE(st.models.empty());
    gemmi::Model &model = st.models[0];

    // Get a chain
    ASSERT_FALSE(model.chains.empty());
    gemmi::Chain &chain = model.chains[0];

    // Get a residue
    ASSERT_FALSE(chain.residues.empty());
    gemmi::Residue &res = chain.residues[0];

    // Get atoms
    ASSERT_FALSE(res.atoms.empty());

    const gemmi::Atom &neighbour_atom = res.atoms[0];

    // Construct atom_overlaps_container_t manually
    coot::atom_overlaps_container_t container;
    container.probe_radius = 1.4; // standard probe radius

    // Populate ligand_atom_neighbour_map with one entry for idx=0
    double vdw_radius = 1.7; // typical carbon vdw radius
    container.ligand_atom_neighbour_map[0].push_back(
        std::make_pair(&neighbour_atom, vdw_radius));

    // Case 1: dot_pt is exactly at the neighbour atom position -> should be inside
    {
        clipper::Coord_orth probe_pos(0, 0, 0);
        clipper::Coord_orth dot_pt(neighbour_atom.pos.x, neighbour_atom.pos.y, neighbour_atom.pos.z);

        bool result = container.is_inside_another_ligand_atom_gemmi(0, probe_pos, dot_pt);
        EXPECT_TRUE(result);
    }

    // Case 2: dot_pt is far away from neighbour atom -> should NOT be inside
    {
        clipper::Coord_orth probe_pos(0, 0, 0);
        clipper::Coord_orth dot_pt(
            neighbour_atom.pos.x + 100,
            neighbour_atom.pos.y + 100,
            neighbour_atom.pos.z + 100);

        bool result = container.is_inside_another_ligand_atom_gemmi(0, probe_pos, dot_pt);
        EXPECT_FALSE(result);
    }

    // Case 3: idx=-1 skips the if block, returns false
    {
        clipper::Coord_orth probe_pos(0, 0, 0);
        clipper::Coord_orth dot_pt(0, 0, 0);

        bool result = container.is_inside_another_ligand_atom_gemmi(-1, probe_pos, dot_pt);
        EXPECT_FALSE(result);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}