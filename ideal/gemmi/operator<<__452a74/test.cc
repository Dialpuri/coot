#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <sstream>
#include "function.hh"

TEST(OracleTest, phi_psi_t) {
    // Case 1: typical phi/psi values (two-arg constructor)
    {
        coot::phi_psi_t pp(120.5f, -45.3f);
        std::ostringstream oss;
        coot::operator_stream_insert_gemmi(oss, pp);
        EXPECT_NEAR(pp.phi, 120.5f, 1e-4);
        EXPECT_NEAR(pp.psi, -45.3f, 1e-4);
        EXPECT_EQ(oss.str(), "120.5 -45.3");
    }

    // Case 2: with tau (three-arg constructor) — tau is stored but NOT printed by operator<<
    {
        coot::phi_psi_t pp(180.0f, 90.0f, 45.0f);
        std::ostringstream oss;
        coot::operator_stream_insert_gemmi(oss, pp);
        EXPECT_NEAR(pp.phi, 180.0f, 1e-4);
        EXPECT_NEAR(pp.psi, 90.0f, 1e-4);
        EXPECT_NEAR(pp.tau, 45.0f, 1e-4);
        EXPECT_EQ(oss.str(), "180 90");
    }

    // Case 3: zero values (edge case)
    {
        coot::phi_psi_t pp(0.0f, 0.0f);
        std::ostringstream oss;
        coot::operator_stream_insert_gemmi(oss, pp);
        EXPECT_NEAR(pp.phi, 0.0f, 1e-4);
        EXPECT_NEAR(pp.psi, 0.0f, 1e-4);
        EXPECT_EQ(oss.str(), "0 0");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}