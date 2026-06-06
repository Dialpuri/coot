#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/mmread.hpp>
#include "api/molecules-container.hh"
#include "function.hh"

// ... TEST(...) blocks here ...

TEST(DifferenceMapPeaksTest, valid_protein_and_map) {
    molecules_container_t mc;

    // Load protein structure
    int imol_protein = mc.read_pdb("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_GE(imol_protein, 0);

    // Load map from MTZ — try DELFWT/PHDELWT (difference map) first
    int imol_map = mc.read_mtz("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.mtz",
                                "DELFWT", "PHDELWT", "", false, false);
    if (imol_map < 0) {
        imol_map = mc.read_mtz("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.mtz",
                                "FWT", "PHWT", "", false, false);
    }
    ASSERT_GE(imol_map, 0);

    float n_rmsd = 3.0f;

    // For the gemmi port, we pass a gemmi::Structure
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Case 1: Valid protein mol, valid map
    {
        std::vector<coot::molecule_t::interesting_place_t> peaks =
            coot::difference_map_peaks_gemmi(mc[imol_map], st, n_rmsd);
        EXPECT_EQ(peaks.size(), 100u);
        EXPECT_FALSE(peaks.empty());
        EXPECT_NEAR(peaks[0].x, 37.6566f, 1e-4);
        EXPECT_NEAR(peaks[0].y, 0.990431f, 1e-4);
        EXPECT_NEAR(peaks[0].z, 56.5793f, 1e-4);
        EXPECT_NEAR(peaks[0].feature_value, -0.867183f, 1e-4);
        EXPECT_NEAR(peaks[0].badness, 40.9118f, 1e-4);
        EXPECT_EQ(peaks[0].button_label, "Index 30 Position: (37.66, 0.99, 56.58) Height -0.87");
    }

    // Case 2: null mol (empty structure) — should produce empty vector with error message
    {
        gemmi::Structure empty_st;
        std::vector<coot::molecule_t::interesting_place_t> peaks =
            coot::difference_map_peaks_gemmi(mc[imol_map], empty_st, n_rmsd);
        EXPECT_EQ(peaks.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}