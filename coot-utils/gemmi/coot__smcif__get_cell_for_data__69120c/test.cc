#include <gtest/gtest.h>
#include <fstream>
#include "function.hh"

TEST(OracleTest, get_cell_for_data) {
    // case: valid CIF file with cell parameters
    {
        std::string test_cif = "/tmp/test-cell.cif";
        std::ofstream ofs(test_cif);
        ofs << "data_test\n";
        ofs << "_cell_length_a  10.0\n";
        ofs << "_cell_length_b  20.0\n";
        ofs << "_cell_length_c  30.0\n";
        ofs << "_cell_angle_alpha  90.0\n";
        ofs << "_cell_angle_beta   90.0\n";
        ofs << "_cell_angle_gamma  90.0\n";
        ofs << "_space_group_IT_number  1\n";
        ofs << "loop_\n";
        ofs << "_atom_site_label\n";
        ofs << "_atom_site_fract_x\n";
        ofs << "_atom_site_fract_y\n";
        ofs << "_atom_site_fract_z\n";
        ofs << "C1 0.5 0.5 0.5\n";
        ofs.close();

        coot::smcif::get_cell_for_data_gemmi(test_cif);
        // Note: This test calls the gemmi port, but the original test
        // used smcif class. The function signature needs to match.
    }

    // case: non-existent file
    {
        std::string cif_file = "/tmp/nonexistent-file.cif";
        coot::smcif::get_cell_for_data_gemmi(cif_file);
        // The function should return a null cell
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}