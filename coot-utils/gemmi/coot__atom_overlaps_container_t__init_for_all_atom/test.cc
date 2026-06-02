#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, atom_overlaps_container_t_gemmi) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    ASSERT_FALSE(st.models.empty());
    gemmi::Model const& model = st.models[0];

    coot::protein_geometry geom;
    geom.init_standard();

    // --- Outputs that mirror atom_overlaps_container_t members ---
    std::vector<gemmi::Residue const*> neighbours;
    std::map<std::string, coot::dictionary_residue_restraints_t> dictionary_map;
    bool have_dictionary = false;
    bool molecule_has_hydrogens = false;
    int overlap_mode = 0;
    unsigned int udd_h_bond_type_handle = 0;
    unsigned int udd_residue_index_handle = 0;

    // Case 1: Construct (equivalent to constructor calling init_for_all_atom internally)
    {
        coot::atom_overlaps_container_t_init_for_all_atom_gemmi(
            neighbours, dictionary_map, have_dictionary, molecule_has_hydrogens, overlap_mode,
            udd_h_bond_type_handle, udd_residue_index_handle,
            model, &geom);

        EXPECT_EQ(overlap_mode, 1);
        EXPECT_EQ(static_cast<int>(neighbours.size()), 268);
        EXPECT_TRUE(have_dictionary);
        EXPECT_TRUE(molecule_has_hydrogens);
        EXPECT_EQ(static_cast<int>(dictionary_map.size()), 21);
        EXPECT_EQ(udd_h_bond_type_handle, 16777219u);
        EXPECT_EQ(udd_residue_index_handle, 16777218u);
    }

    // Case 2: Call init_for_all_atom() explicitly on a fresh set of containers
    // (simulates constructor called with hydrogen_bond=true, then init_for_all_atom())
    {
        std::vector<gemmi::Residue const*> neighbours2;
        std::map<std::string, coot::dictionary_residue_restraints_t> dictionary_map2;
        bool have_dictionary2 = false;
        bool molecule_has_hydrogens2 = false;
        int overlap_mode2 = 0;
        unsigned int udd_h_bond_type_handle2 = 0;
        unsigned int udd_residue_index_handle2 = 0;

        // Simulate constructor for "hydrogen_bond" mode (first init populates neighbours2)
        neighbours2 = neighbours; // pre-populate as constructor would (for "true" case, hydrogen mode)

        int neighbours_before = static_cast<int>(neighbours2.size());

        coot::atom_overlaps_container_t_init_for_all_atom_gemmi(
            neighbours2, dictionary_map2, have_dictionary2, molecule_has_hydrogens2, overlap_mode2,
            udd_h_bond_type_handle2, udd_residue_index_handle2,
            model, &geom);

        int neighbours_after = static_cast<int>(neighbours2.size());

        EXPECT_EQ(neighbours_before, 268);
        EXPECT_EQ(neighbours_after, 536);
        EXPECT_EQ(overlap_mode2, 1);
        EXPECT_TRUE(have_dictionary2);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}