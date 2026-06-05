#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/mmread.hpp>
#include <gemmi/model.hpp>
#include <gemmi/ccp4.hpp>
#include "function.hh"

#include "api/molecules-container.hh"
#include "coot-utils/coot-map-utils.hh"

TEST(OracleTest, mask_map) {
    // case: normal masking with one residue
    {
        molecules_container_t mc;

        // Load the CIF using coot (we still need this for the map)
        int imol = mc.read_pdb("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");
        ASSERT_GE(imol, 0);

        // Load MTZ — this creates a SEPARATE molecule that holds the xmap
        int imol_map = mc.read_mtz("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.mtz", "FWT", "PHWT", "", false, false);
        ASSERT_GE(imol_map, 0);

        coot::molecule_t &mol_map = mc[imol_map];

        // Load the structure with gemmi for residue access
        gemmi::Structure st = gemmi::read_structure_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");

        // Find the LZA residue (chain A, residue 1299)
        const gemmi::Residue *res = nullptr;
        for (gemmi::Model &model : st.models) {
            for (gemmi::Chain &chain : model.chains) {
                if (chain.name == "A") {
                    for (gemmi::Residue &r : chain.residues) {
                        if (r.seqid.num.value == 1299 && r.name == "LZA") {
                            res = &r;
                            break;
                        }
                    }
                }
            }
        }
        ASSERT_NE(res, nullptr);

        // Build neighbour residues vector with just this one residue
        std::vector<const gemmi::Residue *> neighbours;
        neighbours.push_back(res);

        // Get the xmap from the map molecule
        const clipper::Xmap<float> &xmap = mol_map.xmap;

        EXPECT_STREQ(res->name.c_str(), "LZA");
        EXPECT_EQ(res->seqid.num.value, 1299);
        EXPECT_EQ(neighbours.size(), 1u);
        EXPECT_NEAR(xmap.cell().descr().a(), 53.512, 1e-4);

        // Get atom count from the LZA residue to sample
        EXPECT_EQ(res->atoms.size(), 25u);

        // Get first atom position for sampling
        clipper::Coord_orth sample_pt(0, 0, 0);
        if (!res->atoms.empty()) {
            const gemmi::Atom &first_atom = res->atoms[0];
            sample_pt = clipper::Coord_orth(first_atom.pos.x, first_atom.pos.y, first_atom.pos.z);
            std::string padded_name = first_atom.name;
            if (padded_name.size() < 4) padded_name = " " + padded_name;
            EXPECT_STREQ(padded_name.c_str(), " C26");
        }

        // Sample density BEFORE masking
        float density_before = coot::util::density_at_point(xmap, sample_pt);
        EXPECT_NEAR(density_before, 0.754899, 1e-4);

        // Call mask_map
        clipper::Xmap<float> masked = coot::util::mask_map_gemmi(xmap, neighbours);

        // Sample density AFTER masking at the same position
        float density_after = coot::util::density_at_point(masked, sample_pt);
        EXPECT_NEAR(density_after, -10, 1e-4);
        EXPECT_TRUE(density_after != density_before);
        EXPECT_NEAR(masked.cell().descr().a(), 53.512, 1e-4);
    }

    // case: empty neighbours — should return an unmodified copy
    {
        molecules_container_t mc;

        int imol = mc.read_pdb("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");
        ASSERT_GE(imol, 0);
        int imol_map = mc.read_mtz("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.mtz", "FWT", "PHWT", "", false, false);
        ASSERT_GE(imol_map, 0);

        coot::molecule_t &mol_map = mc[imol_map];

        // Load the structure with gemmi for residue access
        gemmi::Structure st = gemmi::read_structure_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");

        // Find the LZA residue (chain A, residue 1299)
        const gemmi::Residue *res = nullptr;
        for (gemmi::Model &model : st.models) {
            for (gemmi::Chain &chain : model.chains) {
                if (chain.name == "A") {
                    for (gemmi::Residue &r : chain.residues) {
                        if (r.seqid.num.value == 1299 && r.name == "LZA") {
                            res = &r;
                            break;
                        }
                    }
                }
            }
        }
        ASSERT_NE(res, nullptr);

        // Get sample point from the residue
        clipper::Coord_orth sample_pt(0, 0, 0);
        if (!res->atoms.empty()) {
            const gemmi::Atom &first_atom = res->atoms[0];
            sample_pt = clipper::Coord_orth(first_atom.pos.x, first_atom.pos.y, first_atom.pos.z);
        }

        const clipper::Xmap<float> &xmap = mol_map.xmap;

        // Empty neighbours
        std::vector<const gemmi::Residue *> empty_neighbours;
        EXPECT_EQ(empty_neighbours.size(), 0u);

        clipper::Xmap<float> unmasked = coot::util::mask_map_gemmi(xmap, empty_neighbours);
        float density_unmasked = coot::util::density_at_point(unmasked, sample_pt);

        float density_before = coot::util::density_at_point(xmap, sample_pt);
        EXPECT_NEAR(density_unmasked, 0.754899, 1e-4);
        EXPECT_TRUE(density_unmasked == density_before);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}