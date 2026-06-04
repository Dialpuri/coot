#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include "function.hh"

TEST(OracleTest, setup_env_residue_atoms_radii_gemmi) {
    // Setup molecules container and load PDB using gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());

    gemmi::Model& model = st.models[0];
    ASSERT_FALSE(model.chains.empty());

    gemmi::Chain& chain = model.chains[0];
    ASSERT_FALSE(chain.residues.empty());

    gemmi::Residue& first_res = chain.residues[0];

    // Verify central residue identity (oracle: A/1 MET)
    EXPECT_EQ(chain.name, "A");
    EXPECT_EQ(first_res.seqid.num.value, 1);
    EXPECT_EQ(first_res.name, "MET");

    // Build neighbour list (first few residues)
    std::vector<gemmi::Residue*> neighbours;
    for (std::size_t i = 0; i < std::min(chain.residues.size(), static_cast<std::size_t>(5)); i++) {
        neighbours.push_back(&chain.residues[i]);
    }
    EXPECT_EQ(neighbours.size(), 5u);

    // Setup protein geometry (default constructor only, matching oracle — no init_standard())
    coot::protein_geometry geom;

    // Setup atom_overlaps_container_t state (using free-function gemmi pattern)
    // Matching the oracle: after init_for_all_atom(), have_dictionary is false
    bool have_dictionary = false;
    std::vector<coot::dictionary_residue_restraints_t> neighb_dictionaries;
    std::map<unsigned long, unsigned int> residue_map;
    std::map<std::string, double> type_to_vdw_radius_map;
    std::vector<double> neighb_atom_radius;

    // Verify oracle outputs
    EXPECT_FALSE(have_dictionary);  // oracle: have_dictionary: false
    EXPECT_EQ(neighb_atom_radius.size(), 0u);  // oracle: neighb_atom_radius_before: 0

    // Count total atoms in model (oracle: 2107)
    std::size_t total_atoms = 0;
    for (auto& ch : model.chains)
        for (auto& res : ch.residues)
            total_atoms += res.atoms.size();

    EXPECT_EQ(total_atoms, 2107u);  // oracle: selection_atom_count: 2107

    // Case 1: Call with all atoms selection
    {
        std::size_t before = neighb_atom_radius.size();
        EXPECT_EQ(before, 0u);  // oracle: neighb_atom_radius_size_before_call: 0

        coot::atom_overlaps_container_t_setup_env_residue_atoms_radii_gemmi(
            neighb_atom_radius,
            have_dictionary,
            neighb_dictionaries,
            residue_map,
            type_to_vdw_radius_map,
            &geom,
            model);

        std::size_t after = neighb_atom_radius.size();
        EXPECT_EQ(after, 2107u);  // oracle: neighb_atom_radius_size_after_call: 2107

        // Count non-zero radii
        int non_zero = 0;
        for (std::size_t i = 0; i < neighb_atom_radius.size(); i++) {
            if (neighb_atom_radius[i] > 0) non_zero++;
        }
        EXPECT_EQ(non_zero, 0);  // oracle: non_zero_radii_count: 0
    }

    // Case 2: Second call — should be no-op due to !empty() guard
    {
        std::size_t before = neighb_atom_radius.size();
        EXPECT_EQ(before, 2107u);  // oracle: second_call_size_before: 2107

        coot::atom_overlaps_container_t_setup_env_residue_atoms_radii_gemmi(
            neighb_atom_radius,
            have_dictionary,
            neighb_dictionaries,
            residue_map,
            type_to_vdw_radius_map,
            &geom,
            model);

        std::size_t after = neighb_atom_radius.size();
        EXPECT_EQ(after, 2107u);  // oracle: second_call_size_after: 2107
        EXPECT_EQ(before, after);  // oracle: second_call_no_op: true
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}