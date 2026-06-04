#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_residue_CA_position_gemmi) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];

    // Case 1: valid residue with CA atom (chain A, residue 10)
    {
        std::string cid = "//A/10/CA";
        std::vector<double> v = coot::get_residue_CA_position_gemmi(cid, model);
        EXPECT_EQ(v.size(), 3u);
        EXPECT_NEAR(v[0], 27.175, 1e-4);
        EXPECT_NEAR(v[1], 12.237, 1e-4);
        EXPECT_NEAR(v[2], 67.71, 1e-4);
    }

    // Case 2: residue that exists but has no CA (e.g. HOH in chain B)
    {
        std::string cid = "//B/1";
        std::vector<double> v = coot::get_residue_CA_position_gemmi(cid, model);
        EXPECT_EQ(v.size(), 0u);
    }

    // Case 3: non-existent residue
    {
        std::string cid = "//A/9999";
        std::vector<double> v = coot::get_residue_CA_position_gemmi(cid, model);
        EXPECT_EQ(v.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}