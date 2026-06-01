#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include <fstream>

TEST(OracleTest, smcif_get_cell) {

    // Case 1: CIF with cell data
    {
        std::string tmp_file = "/tmp/test_cell.cif";
        {
            std::ofstream ofs(tmp_file);
            ofs << "data_test\n";
            ofs << "_cell_length_a  50.0\n";
            ofs << "_cell_length_b  60.0\n";
            ofs << "_cell_length_c  70.0\n";
            ofs << "_cell_angle_alpha  90.0\n";
            ofs << "_cell_angle_beta   90.0\n";
            ofs << "_cell_angle_gamma  120.0\n";
            ofs.close();
        }

        clipper::Cell cell = coot::smcif::get_cell_gemmi(tmp_file);
        EXPECT_NEAR(cell.a(), 50.0, 1e-4);
        EXPECT_NEAR(cell.b(), 60.0, 1e-4);
        EXPECT_NEAR(cell.c(), 70.0, 1e-4);
        EXPECT_NEAR(clipper::Util::rad2d(cell.alpha()), 90.0, 1e-4);
        EXPECT_NEAR(clipper::Util::rad2d(cell.beta()), 90.0, 1e-4);
        EXPECT_NEAR(clipper::Util::rad2d(cell.gamma()), 120.0, 1e-4);
    }

    // Case 2: CIF with no cell data — should throw
    {
        std::string tmp_file = "/tmp/test_no_cell.cif";
        {
            std::ofstream ofs(tmp_file);
            ofs << "data_test\n";
            ofs.close();
        }

        EXPECT_THROW(coot::smcif::get_cell_gemmi(tmp_file), std::runtime_error);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}