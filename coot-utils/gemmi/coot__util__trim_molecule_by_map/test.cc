#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "coot/coot-utils/coot-map-utils.hh"
#include "clipper/core/map_utils.h"

static int count_atoms_in_model(const gemmi::Model &model) {
    int count = 0;
    for (const gemmi::Chain &chain : model.chains) {
        for (const gemmi::Residue &res : chain.residues) {
            count += static_cast<int>(res.atoms.size());
        }
    }
    return count;
}

TEST(OracleTest, trim_molecule_by_map) {
    // Load PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());
    gemmi::Model &model = st.models[0];

    // Create map from MTZ
    clipper::Xmap<float> xmap;
    bool map_ok = coot::util::map_fill_from_mtz(&xmap,
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.mtz",
        "FWT", "PHWT", "", 0);
    ASSERT_TRUE(map_ok);

    // Compute map level
    clipper::Map_stats stats(xmap);
    float sigma_level = 3.0f;
    float map_level = stats.mean() + sigma_level * stats.std_dev();

    // Case 1: waters only, DELETE mode
    {
        int atoms_before = count_atoms_in_model(model);

        short int remove_flag = coot::util::TRIM_BY_MAP_DELETE;
        short int waters_only = 1;

        int n_changed = coot::util::trim_molecule_by_map_gemmi(model, xmap, map_level, remove_flag, waters_only);

        int atoms_after = count_atoms_in_model(model);

        EXPECT_NEAR(map_level, 0.904088f, 1e-4);
        EXPECT_EQ(atoms_before, 2107);
        EXPECT_EQ(atoms_after, 2106);
        EXPECT_EQ(n_changed, 1);
    }

    // Case 2: non-waters only (waters_only=0), ZERO_OCC mode, very high threshold
    {
        int atoms_before = count_atoms_in_model(model);

        short int remove_flag = coot::util::TRIM_BY_MAP_ZERO_OCC;
        short int waters_only = 0;

        // Use a very high map level so nothing gets trimmed
        float high_level = stats.mean() + 10.0f * stats.std_dev();
        int n_changed = coot::util::trim_molecule_by_map_gemmi(model, xmap, high_level, remove_flag, waters_only);

        int atoms_after = count_atoms_in_model(model);

        EXPECT_NEAR(high_level, 3.01363f, 1e-4);
        EXPECT_EQ(atoms_before, 2106);
        EXPECT_EQ(atoms_after, 2106);
        EXPECT_EQ(n_changed, 2106);
    }

    // Case 3: complementary - DELETE mode, non-waters (waters_only=0), high threshold to actually delete
    {
        // Reload fresh molecule for this case to get a clean state
        gemmi::Structure st2 = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        ASSERT_FALSE(st2.models.empty());
        gemmi::Model &model2 = st2.models[0];

        int atoms_before = count_atoms_in_model(model2);

        short int remove_flag = coot::util::TRIM_BY_MAP_DELETE;
        short int waters_only = 0;

        // Use a very high map level so atoms actually get deleted
        float high_level = stats.mean() + 10.0f * stats.std_dev();
        int n_changed = coot::util::trim_molecule_by_map_gemmi(model2, xmap, high_level, remove_flag, waters_only);

        int atoms_after = count_atoms_in_model(model2);

        EXPECT_EQ(atoms_before, 2107);
        EXPECT_TRUE(atoms_after < atoms_before) << "DELETE mode should reduce atom count: before=" << atoms_before << " after=" << atoms_after;
        EXPECT_EQ(n_changed, atoms_before - atoms_after);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
