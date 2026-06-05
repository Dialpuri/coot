#include <gtest/gtest.h>
#include <gemmi/model.hpp>
#include <gemmi/cif.hpp>
#include <gemmi/mmcif.hpp>
#include <clipper/core/xmap.h>
#include <clipper/core/cell.h>
#include <clipper/core/spacegroup.h>
#include <clipper/core/coords.h>
#include "function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"
#include <vector>
#include <cmath>
#include <unistd.h>

TEST(OracleTest, make_map_mask_gemmi) {
    gemmi::cif::Document doc = gemmi::cif::read_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");
    gemmi::Structure st = gemmi::make_structure(std::move(doc));
    ASSERT_FALSE(st.models.empty());

    // Construct clipper types from gemmi structure
    clipper::Spgr_descr spgr_descr(st.spacegroup_hm.c_str());
    clipper::Spacegroup spacegroup(spgr_descr);
    clipper::Cell_descr cell_descr(st.cell.a, st.cell.b, st.cell.c,
                                   st.cell.alpha, st.cell.beta, st.cell.gamma);
    clipper::Cell cell(cell_descr);

    // Grid dimensions: nu=108, nv=144, nw=144 (from example.mtz)
    clipper::Grid_sampling gs(108, 144, 144);

    // Verify setup values
    EXPECT_EQ(spacegroup.symbol_hm(), "P 21 21 21");
    EXPECT_NEAR(cell.a(), 53.512, 1e-4);
    EXPECT_NEAR(cell.b(), 71.311, 1e-4);
    EXPECT_NEAR(cell.c(), 72.101, 1e-4);
    EXPECT_EQ(gs.nu(), 108);
    EXPECT_EQ(gs.nv(), 144);
    EXPECT_EQ(gs.nw(), 144);

    // Select atoms from chain A, residue 1299 (LZA ligand)
    gemmi::Model& model = st.models[0];
    std::vector<const gemmi::Atom*> selected_atoms;
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == "A") {
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == 1299) {
                    for (gemmi::Atom& atom : res.atoms) {
                        selected_atoms.push_back(&atom);
                    }
                }
            }
        }
    }

    // Case 1: smooth = 0.0f (no smoothing), selecting LZA ligand atoms
    {
        EXPECT_EQ(static_cast<int>(selected_atoms.size()), 25);

        float radius = 3.0f;
        float smooth = 0.0f;

        clipper::Xmap<float> result = coot::util::make_map_mask_gemmi(
            spacegroup, cell, gs, selected_atoms, radius, smooth);

        clipper::Grid_sampling res_gs = result.grid_sampling();
        EXPECT_EQ(res_gs.nu(), 108);
        EXPECT_EQ(res_gs.nv(), 144);
        EXPECT_EQ(res_gs.nw(), 144);

        int count_one = 0;
        int count_outer = 0;
        int count_other = 0;
        float max_val = -999;
        float min_val = 999;
        for (clipper::Xmap_base::Map_reference_index ix = result.first(); !ix.last(); ix.next()) {
            float v = result[ix];
            if (v == 1.0f) count_one++;
            else if (v < -1.0f) count_outer++;
            else count_other++;
            if (v > max_val) max_val = v;
            if (v < min_val) min_val = v;
        }
        EXPECT_EQ(count_one, 6777);
        EXPECT_EQ(count_outer, 553095);
        EXPECT_EQ(count_other, 0);
        EXPECT_NEAR(max_val, 1.0, 1e-4);
        EXPECT_NEAR(min_val, -1.1, 1e-4);
    }

    // Case 2: smooth = 1.0f (with smoothing)
    {
        EXPECT_EQ(static_cast<int>(selected_atoms.size()), 25);

        float radius = 3.0f;
        float smooth = 1.0f;

        clipper::Xmap<float> result = coot::util::make_map_mask_gemmi(
            spacegroup, cell, gs, selected_atoms, radius, smooth);

        int count_one = 0;
        int count_outer = 0;
        int count_other = 0;
        float max_val = -999;
        float min_val = 999;
        for (clipper::Xmap_base::Map_reference_index ix = result.first(); !ix.last(); ix.next()) {
            float v = result[ix];
            if (v == 1.0f) count_one++;
            else if (v < -1.0f) count_outer++;
            else count_other++;
            if (v > max_val) max_val = v;
            if (v < min_val) min_val = v;
        }
        EXPECT_EQ(count_one, 2032);
        EXPECT_EQ(count_outer, 554182);
        EXPECT_EQ(count_other, 3658);
        EXPECT_NEAR(max_val, 1.0, 1e-4);
        EXPECT_NEAR(min_val, -1.1, 1e-4);
    }

    // Case 3: empty selection (non-existent residue) — should produce all-outer map
    {
        std::vector<const gemmi::Atom*> empty_atoms;
        EXPECT_EQ(static_cast<int>(empty_atoms.size()), 0);

        float radius = 3.0f;
        float smooth = 0.0f;

        clipper::Xmap<float> result = coot::util::make_map_mask_gemmi(
            spacegroup, cell, gs, empty_atoms, radius, smooth);

        int count_one = 0;
        int count_outer = 0;
        for (clipper::Xmap_base::Map_reference_index ix = result.first(); !ix.last(); ix.next()) {
            float v = result[ix];
            if (v == 1.0f) count_one++;
            else if (v < -1.0f) count_outer++;
        }
        EXPECT_EQ(count_one, 0);
        EXPECT_EQ(count_outer, 559872);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}