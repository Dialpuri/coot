#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// case: oracle test for standard_deviation_temperature_factor
TEST(OracleTest, StandardDeviationTemperatureFactor) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_GT(st.models.size(), 0) << "PDB loading failed";

    // case 1: valid molecule with all atoms (no cutoffs)
    {
        float sd = coot::util::standard_deviation_temperature_factor_gemmi(
            st.models[0], 0.0f, 999.0f, false, false);
        EXPECT_NEAR(sd, 15.3039, 1e-4);
    }

    // case 2: apply low cutoff (exclude low B-factors)
    {
        float sd = coot::util::standard_deviation_temperature_factor_gemmi(
            st.models[0], 50.0f, 999.0f, true, false);
        EXPECT_NEAR(sd, 22.6694, 1e-4);
    }

    // case 3: apply high cutoff (exclude high B-factors)
    {
        float sd = coot::util::standard_deviation_temperature_factor_gemmi(
            st.models[0], 0.0f, 30.0f, false, true);
        EXPECT_NEAR(sd, 25.0553, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}