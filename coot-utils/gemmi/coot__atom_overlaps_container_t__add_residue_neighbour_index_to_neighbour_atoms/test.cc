#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, add_residue_neighbour_index_to_neighbour_atoms) {
    // Load PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    ASSERT_GE(static_cast<int>(st.models.size()), 1);

    // Get two residues from chain A to use as neighbours
    gemmi::Model* model = &st.models[0];
    ASSERT_NE(model, nullptr);

    // Find chain A
    gemmi::Chain* chain = model->find_chain("A");
    ASSERT_NE(chain, nullptr);

    // Pick residue index 10 and 11 from chain A (MMDB GetResidue is 0-based)
    ASSERT_LT(10u, chain->residues.size());
    ASSERT_LT(11u, chain->residues.size());
    gemmi::Residue* res10 = &chain->residues[10];
    gemmi::Residue* res11 = &chain->residues[11];

    ASSERT_NE(res10, nullptr);
    ASSERT_NE(res11, nullptr);

    // Create and configure the container
    coot::atom_overlaps_container_t container;
    container.neighbours.push_back(res10);
    container.neighbours.push_back(res11);

    // (inputs verified by successful construction above)

    // Before call: check udd handle is unset
    EXPECT_EQ(container.udd_residue_index_handle, 0);

    // Call the function
    container.add_residue_neighbour_index_to_neighbour_atoms_gemmi();

    // After call: check udd handle
    EXPECT_EQ(container.udd_residue_index_handle, 16777217);

    // Verify atoms in neighbour residues have the correct UDD data
    int natoms_10 = static_cast<int>(res10->atoms.size());
    int natoms_11 = static_cast<int>(res11->atoms.size());

    EXPECT_EQ(natoms_10, 4);
    EXPECT_EQ(natoms_11, 9);

    // Check a sample atom from residue 0 has UDD value 0
    int sample_udd_0 = -999;
    if (natoms_10 > 0) {
        sample_udd_0 = container.atom_neighbour_index[&res10->atoms[0]];
    }
    EXPECT_EQ(sample_udd_0, 0);

    // Check a sample atom from residue 1 has UDD value 1
    int sample_udd_1 = -999;
    if (natoms_11 > 0) {
        sample_udd_1 = container.atom_neighbour_index[&res11->atoms[0]];
    }
    EXPECT_EQ(sample_udd_1, 1);

    // Edge case: empty neighbours vector
    {
        coot::atom_overlaps_container_t container2;
        // neighbours is empty
        EXPECT_EQ(container2.neighbours.size(), 0u);
        container2.add_residue_neighbour_index_to_neighbour_atoms_gemmi();
        EXPECT_EQ(container2.udd_residue_index_handle, 16777217);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}