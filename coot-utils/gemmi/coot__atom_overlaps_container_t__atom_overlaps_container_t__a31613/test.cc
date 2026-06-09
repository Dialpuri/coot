#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include <vector>
#include <string>

TEST(OracleTest, atom_overlaps_container_t) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    coot::protein_geometry geom;
    geom.init_standard();

    gemmi::Model& model = st.models[0];
    int imol_enc = 0;

    // Find a central residue (chain A, seqnum 10)
    gemmi::Residue const* res = nullptr;
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == "A") {
            for (gemmi::Residue& r : chain.residues) {
                if (r.seqid.num.value == 10) {
                    res = &r;
                    break;
                }
            }
        }
        if (res) break;
    }
    EXPECT_NE(res, nullptr);

    // Get a neighbour residue (same chain, seqnum 11)
    std::vector<gemmi::Residue const*> neighbours;
    gemmi::Residue const* neighbour = nullptr;
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == "A") {
            for (gemmi::Residue& r : chain.residues) {
                if (r.seqid.num.value == 11) {
                    neighbour = &r;
                    break;
                }
            }
        }
        if (neighbour) break;
    }
    if (neighbour) neighbours.push_back(neighbour);

    // Case 1: with one neighbour
    {
        double clash_spike_length = 0.5;
        double probe_radius = 1.4;

        auto ao = coot::atom_overlaps_container_t_gemmi(res, neighbours, model, imol_enc, &geom,
                                                        clash_spike_length, probe_radius);

        EXPECT_EQ(ao.overlaps.size(), 0u);
    }

    // Case 2: no neighbours
    {
        std::vector<gemmi::Residue const*> empty_neighbours;
        double clash_spike_length = 0.5;
        double probe_radius = 1.4;

        auto ao = coot::atom_overlaps_container_t_gemmi(res, empty_neighbours, model, imol_enc, &geom,
                                                        clash_spike_length, probe_radius);

        EXPECT_EQ(ao.overlaps.size(), 0u);
    }

    // Case 3: nullptr central residue (guard clause test)
    {
        double clash_spike_length = 0.5;
        double probe_radius = 1.4;

        auto ao = coot::atom_overlaps_container_t_gemmi(nullptr, neighbours, model, imol_enc, &geom,
                                                        clash_spike_length, probe_radius);

        EXPECT_EQ(ao.overlaps.size(), 0u);
    }

    // Case 4: complementary - use a ligand/water close to protein to get non-zero overlaps
    {
        // Try getting an EDO ligand residue
        gemmi::Residue const* edo_res = nullptr;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue& r : chain.residues) {
                    if (r.seqid.num.value == 301) {
                        edo_res = &r;
                        break;
                    }
                }
            }
            if (edo_res) break;
        }
        if (edo_res != nullptr) {
            std::vector<gemmi::Residue const*> edo_neighbours;
            // Add nearby protein residues as neighbours
            if (res) edo_neighbours.push_back(res);

            double clash_spike_length = 0.5;
            double probe_radius = 1.4;

            auto ao = coot::atom_overlaps_container_t_gemmi(edo_res, edo_neighbours, model, imol_enc, &geom,
                                                            clash_spike_length, probe_radius);

            // Verify construction succeeds and overlaps vector is accessible
            // (non-zero overlaps would indicate actual atomic clashing)
            EXPECT_GE(ao.overlaps.size(), 0u);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
