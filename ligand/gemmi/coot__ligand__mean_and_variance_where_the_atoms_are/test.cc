#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "coot/api/molecules-container.hh"
#include "function.hh"

TEST(OracleTest, mean_and_variance_where_the_atoms_are_gemmi) {
    molecules_container_t mc;
    mc.geometry_init_standard();

    int imol_pdb = mc.read_pdb("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_GE(imol_pdb, 0);

    int imol_map = mc.read_mtz("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.mtz",
                                "FWT", "PHWT", "", false, false);
    ASSERT_GE(imol_map, 0);

    const clipper::Xmap<float> &xmap = mc[imol_map].xmap;

    // Load the same PDB with gemmi to get gemmi::Model
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());

    const gemmi::Model &model = st.models[0];

    // Case 1: valid molecule + valid map
    {
        std::pair<float, float> mv = coot::ligand::mean_and_variance_where_the_atoms_are_gemmi(model, xmap);

        EXPECT_FALSE(std::isnan(mv.first));
        EXPECT_FALSE(std::isinf(mv.first));
        EXPECT_GT(mv.first, 0.0f);
        EXPECT_LT(mv.first, 2.0f);

        EXPECT_FALSE(std::isnan(mv.second));
        EXPECT_FALSE(std::isinf(mv.second));
        EXPECT_GT(mv.second, 0.0f);
        EXPECT_LT(mv.second, 2.0f);
    }

    // Case 2: same molecule + map, different call
    {
        std::pair<float, float> mv = coot::ligand::mean_and_variance_where_the_atoms_are_gemmi(model, xmap);

        EXPECT_FALSE(std::isnan(mv.first));
        EXPECT_FALSE(std::isinf(mv.first));
        EXPECT_GT(mv.first, 0.0f);
        EXPECT_LT(mv.first, 2.0f);

        EXPECT_FALSE(std::isnan(mv.second));
        EXPECT_FALSE(std::isinf(mv.second));
        EXPECT_GT(mv.second, 0.0f);
        EXPECT_LT(mv.second, 2.0f);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}