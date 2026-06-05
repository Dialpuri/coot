#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include <clipper/clipper.h>
#include <vector>

TEST(OracleTestGemmi, fffear_search_fill_nxmap) {
    // Load PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Create a selection of all atoms in chain A
    std::vector<const gemmi::Atom*> atom_selection;
    gemmi::Model& model = st.models[0];
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == "A") {
            for (gemmi::Residue& res : chain.residues) {
                for (gemmi::Atom& atom : res.atoms) {
                    atom_selection.push_back(&atom);
                }
            }
        }
    }

    int n_atoms = static_cast<int>(atom_selection.size());

    ASSERT_GT(n_atoms, 0);
    EXPECT_EQ(n_atoms, 2106);

    // Compute mid_point of the selection
    float sum_x = 0, sum_y = 0, sum_z = 0;
    for (int i = 0; i < n_atoms; i++) {
        if (atom_selection[i]) {
            sum_x += static_cast<float>(atom_selection[i]->pos.x);
            sum_y += static_cast<float>(atom_selection[i]->pos.y);
            sum_z += static_cast<float>(atom_selection[i]->pos.z);
        }
    }
    float mid_x = sum_x / n_atoms;
    float mid_y = sum_y / n_atoms;
    float mid_z = sum_z / n_atoms;
    clipper::Coord_orth mid_point(mid_x, mid_y, mid_z);

    // Create an Xmap to satisfy the fffear_search constructor
    clipper::Cell_descr cell_descr(80.0, 80.0, 80.0, 90.0, 90.0, 90.0);
    clipper::Cell cell(cell_descr);
    clipper::Spacegroup sg = clipper::Spacegroup::p1();
    clipper::Resolution reso(2.0);
    clipper::Grid_sampling grid_sampling(sg, cell, reso);
    clipper::Grid_range grid_extent(cell, grid_sampling, 50.0);

    clipper::Xmap<float> xmap;
    xmap.init(sg, cell, grid_sampling);
    // Set xmap to 1.0 everywhere so it's not empty
    clipper::Xmap<float>::Map_reference_index ix;
    for (ix = xmap.first(); !ix.last(); ix.next()) {
        xmap[ix] = 1.0f;
    }

    // Create NXmap from xmap (mimicking fffear_search constructor)
    clipper::NXmap<float> nxmap;
    nxmap.init(xmap.cell(), grid_sampling, grid_extent);

    // First reset nxmap to zeros.
    nxmap = 0.0f;

    int n_points_before = 0;
    clipper::NXmap<float>::Map_reference_index inx;
    for (inx = nxmap.first(); !inx.last(); inx.next()) {
        if (nxmap[inx] != 0.0f) n_points_before++;
    }

    // Call fill_nxmap_gemmi
    int n_points = coot::util::fffear_search::fill_nxmap_gemmi(nxmap, atom_selection, mid_point);

    int n_points_after = 0;
    for (inx = nxmap.first(); !inx.last(); inx.next()) {
        if (nxmap[inx] != 0.0f) n_points_after++;
    }

    EXPECT_EQ(n_points, 1535274);
    EXPECT_EQ(n_points_before, 0);
    EXPECT_EQ(n_points_after, 210108);

    // Case 2: Empty selection — should return 0 points
    {
        std::vector<const gemmi::Atom*> atom_selection2;
        // Try to find residue with seqid 99999 in chain A (won't exist)
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue& res : chain.residues) {
                    if (res.seqid.num.value == 99999) {
                        for (gemmi::Atom& atom : res.atoms) {
                            atom_selection2.push_back(&atom);
                        }
                    }
                }
            }
        }

        int n_atoms2 = static_cast<int>(atom_selection2.size());

        EXPECT_EQ(n_atoms2, 0);

        nxmap = 0.0f;
        int n_points2 = coot::util::fffear_search::fill_nxmap_gemmi(nxmap, atom_selection2, mid_point);

        int n_points2_after = 0;
        for (inx = nxmap.first(); !inx.last(); inx.next()) {
            if (nxmap[inx] != 0.0f) n_points2_after++;
        }

        EXPECT_EQ(n_points2, 0);
        EXPECT_EQ(n_points2_after, 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}