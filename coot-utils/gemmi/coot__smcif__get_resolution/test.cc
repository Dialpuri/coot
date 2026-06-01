#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include <iostream>
#include <fstream>
#include <cmath>

TEST(OracleTest, get_resolution) {
    // Construct a clipper::Cell with typical dimensions
    clipper::Cell cell(clipper::Cell_descr(50.0, 60.0, 70.0, 90.0, 90.0, 90.0));

    // Case 1: Valid CIF with reflection data
    {
        std::string tmp_file = "/tmp/test_refln.cif";
        std::ofstream ofs(tmp_file);
        ofs << "data_test\n"
            << "loop_\n"
            << "_refln_index_h\n"
            << "_refln_index_k\n"
            << "_refln_index_l\n"
            << "1 0 0\n"
            << "0 2 0\n"
            << "0 0 3\n"
            << "\n";
        ofs.close();

        clipper::Resolution reso = coot::smcif::get_resolution_gemmi(cell, tmp_file);

        EXPECT_FALSE(reso.is_null());
        EXPECT_NEAR(reso.limit(), 23.3333, 1e-4);
    }

    // Case 2: Non-existent file (should still return non-null Resolution with inf limit)
    {
        std::string file_name = "/nonexistent/file.cif";
        clipper::Resolution reso = coot::smcif::get_resolution_gemmi(cell, file_name);

        EXPECT_FALSE(reso.is_null());
        EXPECT_TRUE(std::isinf(reso.limit()));
    }

    // Case 3: CIF without reflection data
    {
        std::string file_name = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif";
        clipper::Resolution reso = coot::smcif::get_resolution_gemmi(cell, file_name);

        EXPECT_FALSE(reso.is_null());
        EXPECT_TRUE(std::isinf(reso.limit()));
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}