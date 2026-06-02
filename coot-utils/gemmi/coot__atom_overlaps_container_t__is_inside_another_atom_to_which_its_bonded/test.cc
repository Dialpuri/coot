#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, is_inside_another_atom_to_which_its_bonded) {
    // Load a PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Collect all atoms from the first model
    gemmi::Model &model = st.models[0];
    std::vector<const gemmi::Atom *> all_atoms;
    for (gemmi::Chain &chain : model.chains) {
        for (gemmi::Residue &residue : chain.residues) {
            for (gemmi::Atom &atom : residue.atoms) {
                all_atoms.push_back(&atom);
            }
        }
    }

    ASSERT_GE(static_cast<int>(all_atoms.size()), 2);

    // Get first two atoms from the same residue (they are likely bonded, e.g. N-CA)
    const gemmi::Atom *at_center = all_atoms[0];
    const gemmi::Atom *at_neighbor = all_atoms[1];

    // Build atom_selection array (size 2: index 0 = center, index 1 = neighbor)
    const gemmi::Atom **atom_selection = new const gemmi::Atom *[2];
    atom_selection[0] = at_center;
    atom_selection[1] = at_neighbor;

    // VDW radii (use typical values: C=1.7, N=1.55, O=1.52, S=1.8)
    std::vector<double> neighb_atom_radius(2, 1.7);
    {
        std::string elem = at_center->element.name();
        if (elem == "N") neighb_atom_radius[0] = 1.55;
        else if (elem == "O") neighb_atom_radius[0] = 1.52;
        else if (elem == "S") neighb_atom_radius[0] = 1.8;
    }
    {
        std::string elem = at_neighbor->element.name();
        if (elem == "N") neighb_atom_radius[1] = 1.55;
        else if (elem == "O") neighb_atom_radius[1] = 1.52;
        else if (elem == "S") neighb_atom_radius[1] = 1.8;
    }

    // Bonded neighbor indices: atom 0 has neighbor at index 1
    std::vector<int> bonded_indices = {1};

    // Create dummy container to call the method
    coot::atom_overlaps_container_t container;

    // === CASE 1: Surface point IS inside bonded neighbor (should return true) ===
    // Place surface point AT the neighbor's position (distance = 0 < r_2^2)
    clipper::Coord_orth neighbor_pos = coot::co_gemmi(at_neighbor);
    {
        int atom_idx = 0;
        bool result = container.is_inside_another_atom_to_which_its_bonded_gemmi(
            atom_idx, at_center, neighbor_pos, bonded_indices, atom_selection, neighb_atom_radius);
        EXPECT_TRUE(result);
    }

    // === CASE 2: Surface point is FAR from bonded neighbor (should return false) ===
    // Place surface point 100 Angstroms away
    clipper::Coord_orth far_point(neighbor_pos.x() + 100.0, neighbor_pos.y(), neighbor_pos.z());
    {
        int atom_idx = 0;
        bool result = container.is_inside_another_atom_to_which_its_bonded_gemmi(
            atom_idx, at_center, far_point, bonded_indices, atom_selection, neighb_atom_radius);
        EXPECT_FALSE(result);
    }

    // === CASE 3: Empty bonded list (should return false) ===
    std::vector<int> empty_bonded;
    {
        int atom_idx = 0;
        bool result = container.is_inside_another_atom_to_which_its_bonded_gemmi(
            atom_idx, at_center, neighbor_pos, empty_bonded, atom_selection, neighb_atom_radius);
        EXPECT_FALSE(result);
    }

    delete[] atom_selection;
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}