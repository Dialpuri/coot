#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/calculate.hpp>
#include <filesystem>
#include <string>
#include "function.hh"

TEST(OracleTest, write_coordinates_gemmi) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    int n_selected = static_cast<int>(gemmi::count_atom_sites(st));
    EXPECT_EQ(n_selected, 2107);

    // Case 1: Write as PDB (no extension match for .cif or shelx)
    {
        std::string file_name = "/tmp/oracle_write_coords_test.pdb";
        int err = coot::write_coordinates_gemmi(st, file_name, n_selected);
        EXPECT_EQ(err, 0);
        bool file_exists = std::filesystem::exists(file_name);
        EXPECT_TRUE(file_exists);
        if (file_exists) {
            std::uintmax_t sz = std::filesystem::file_size(file_name);
            EXPECT_GT(sz, 0u);
            EXPECT_EQ(sz, 175381);
        }
    }

    // Case 2: Write as CIF
    {
        std::string file_name = "/tmp/oracle_write_coords_test.cif";
        int err = coot::write_coordinates_gemmi(st, file_name, n_selected);
        EXPECT_EQ(err, 0);
        bool file_exists = std::filesystem::exists(file_name);
        EXPECT_TRUE(file_exists);
        if (file_exists) {
            std::uintmax_t sz = std::filesystem::file_size(file_name);
            EXPECT_GT(sz, 0u);
            EXPECT_EQ(sz, 146515);
        }
    }

    // Case 3: Empty selection — n_selected_atoms == 0 should guard and return 1
    {
        std::string file_name = "/tmp/oracle_write_coords_empty.pdb";
        int err = coot::write_coordinates_gemmi(st, file_name, 0);
        EXPECT_EQ(err, 1);
        bool file_exists = std::filesystem::exists(file_name);
        EXPECT_FALSE(file_exists);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}