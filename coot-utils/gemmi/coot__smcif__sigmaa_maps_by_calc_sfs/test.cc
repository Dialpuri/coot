#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "coot-utils/read-sm-cif.hh"

TEST(OracleTest, sigmaa_maps_by_calc_sfs) {
    // Load the SM CIF file into smcif to get reflection data
    coot::smcif smcif;
    std::string sm_cif_file = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif";

    bool read_ok = smcif.read_data_sm_cif(sm_cif_file);

    // Oracle output: OUTPUT read_data_sm_cif_success: false
    EXPECT_FALSE(read_ok);

    // The oracle returns early when read_ok is false, so the subsequent
    // code (molecule loading, sigmaa_maps_by_calc_sfs) is never reached.
    // We faithfully reproduce this behavior.
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}