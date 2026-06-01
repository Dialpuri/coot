#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/cifdoc.hpp>
#include <gemmi/read_cif.hpp>
#include "function.hh"
#include <fstream>

TEST(OracleTest, get_cell_for_data) {
    // Create a temporary CIF file with full cell data
    const char* tmp_cif = "/tmp/test_cell.cif";
    {
        std::ofstream ofs(tmp_cif);
        ofs << "data_test\n";
        ofs << "_cell_length_a  10.0\n";
        ofs << "_cell_length_b  20.0\n";
        ofs << "_cell_length_c  30.0\n";
        ofs << "_cell_angle_alpha  90.0\n";
        ofs << "_cell_angle_beta   90.0\n";
        ofs << "_cell_angle_gamma  90.0\n";
        ofs.close();
    }

    // Case 1: Valid CIF with cell data
    {
        gemmi::cif::Document doc = gemmi::read_cif_gz(tmp_cif);

        EXPECT_FALSE(doc.blocks.empty()); // read_ok: true

        const gemmi::cif::Block& data = doc.blocks[0];
        clipper::Cell cell = coot::smcif::get_cell_for_data_gemmi(data);
        EXPECT_FALSE(cell.is_null()); // cell_is_null: false

        // cell_a: 10, cell_b: 20, cell_c: 30
        EXPECT_NEAR(cell.a(), 10.0, 1e-4);
        EXPECT_NEAR(cell.b(), 20.0, 1e-4);
        EXPECT_NEAR(cell.c(), 30.0, 1e-4);

        // Angles: oracle printed in degrees
        double deg = 180.0 / 3.14159265358979323846;
        EXPECT_NEAR(cell.alpha() * deg, 90.0, 1e-4); // cell_alpha: 90
        EXPECT_NEAR(cell.beta() * deg, 90.0, 1e-4);   // cell_beta: 90
        EXPECT_NEAR(cell.gamma() * deg, 90.0, 1e-4);  // cell_gamma: 90
    }

    // Case 2: Empty data object (no cell parameters loaded)
    {
        gemmi::cif::Block data;

        clipper::Cell cell = coot::smcif::get_cell_for_data_gemmi(data);
        EXPECT_TRUE(cell.is_null()); // cell_is_null: true
    }

    // Case 3: CIF with missing cell parameters (partial data)
    {
        const char* tmp_cif2 = "/tmp/test_cell_partial.cif";
        {
            std::ofstream ofs(tmp_cif2);
            ofs << "data_test_partial\n";
            ofs << "_cell_length_a  15.0\n";
            ofs << "_cell_length_b  25.0\n";
            // missing c, angles
            ofs.close();
        }

        gemmi::cif::Document doc = gemmi::read_cif_gz(tmp_cif2);

        EXPECT_FALSE(doc.blocks.empty()); // read_ok: true

        const gemmi::cif::Block& data = doc.blocks[0];
        clipper::Cell cell = coot::smcif::get_cell_for_data_gemmi(data);
        EXPECT_TRUE(cell.is_null()); // cell_is_null: true
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}