#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include "geometry/protein-geometry.hh"

TEST(OracleTest, atom_overlaps_container_t) {
    coot::protein_geometry geom;
    geom.init_standard();

    // Load the PDB with gemmi instead of molecules_container_t
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_EQ(st.models.size(), 1u);
    gemmi::Model const& model = st.models[0];

    // Case 1: standard construction with default probe_radius and clash_spike_length
    {
        bool ignore_water_contacts_flag = false;
        double clash_spike_length_in = 0.5;
        double probe_radius_in = 1.4;

        auto ao = coot::atom_overlaps_container_t_gemmi(
            model, &geom, ignore_water_contacts_flag, clash_spike_length_in, probe_radius_in);

        EXPECT_EQ(ao.neighbours.size(), 268u);
        EXPECT_TRUE(ao.have_dictionary);
        EXPECT_FALSE(ao.molecule_has_hydrogens);
        EXPECT_EQ(ao.overlap_mode, 1);
    }

    // Case 2: with ignore_water_contacts_flag = true
    {
        bool ignore_water_contacts_flag = true;
        double clash_spike_length_in = 0.3;
        double probe_radius_in = 1.0;

        auto ao = coot::atom_overlaps_container_t_gemmi(
            model, &geom, ignore_water_contacts_flag, clash_spike_length_in, probe_radius_in);

        EXPECT_EQ(ao.neighbours.size(), 268u);
        EXPECT_TRUE(ao.have_dictionary);
        EXPECT_TRUE(ao.ignore_water_contacts_flag);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}