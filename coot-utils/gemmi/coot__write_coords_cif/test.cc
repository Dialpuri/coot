#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/mmread.hpp>
#include <gemmi/calculate.hpp>
#include "function.hh"
#include <fstream>
#include <iostream>
#include <string>

static bool file_exists_and_size(const std::string &path, long &sz) {
    std::ifstream f(path.c_str(), std::ios::binary | std::ios::ate);
    if (!f.is_open()) { sz = 0; return false; }
    sz = (long)f.tellg();
    f.close();
    return true;
}

TEST(OracleTest, write_coords_cif) {
    gemmi::Structure st = gemmi::read_structure_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");

    size_t n_atoms = gemmi::count_atom_sites(st);
    EXPECT_EQ(n_atoms, 2474);

    // case: valid molecule write
    {
        std::string out_file = "/tmp/oracle_test_write_cif_1.cif";
        int result = coot::write_coords_cif_gemmi(st, out_file);

        long file_size = 0;
        bool file_exists = file_exists_and_size(out_file, file_size);

        EXPECT_EQ(gemmi::count_atom_sites(st), 2474);
        EXPECT_EQ(result, 0);
        EXPECT_TRUE(file_exists);
        EXPECT_EQ(file_size, 181286L);
    }

    // case: write to a different path
    {
        std::string out_file = "/tmp/oracle_test_write_cif_2.cif";
        int result = coot::write_coords_cif_gemmi(st, out_file);

        long file_size = 0;
        bool file_exists = file_exists_and_size(out_file, file_size);

        EXPECT_EQ(gemmi::count_atom_sites(st), 2474);
        EXPECT_EQ(result, 0);
        EXPECT_TRUE(file_exists);
        EXPECT_EQ(file_size, 181286L);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}