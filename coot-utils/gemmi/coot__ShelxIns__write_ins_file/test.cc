#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include <fstream>
#include <string>

TEST(OracleTest, write_ins_file) {
    // Load PDB into gemmi::Structure
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Set cell and space group (PDB file has no crystallographic info)
    st.cell = gemmi::UnitCell(50.0, 60.0, 70.0, 90.0, 90.0, 90.0);
    st.spacegroup_hm = "P 1";

    // Case 1: mol_is_from_shelx_ins = false
    {
        std::string outfile = "/tmp/test_shelx_ins_1.ins";
        coot::ShelxIns shelx;
        std::pair<int, std::string> result = shelx.write_ins_file_gemmi(st, outfile, false);

        EXPECT_EQ(result.first, 1);
        EXPECT_EQ(result.second, "INFO:: SHELXL file /tmp/test_shelx_ins_1.ins written.");

        std::ifstream check(outfile);
        bool file_exists = check.good();
        EXPECT_TRUE(file_exists);

        std::string first_line;
        if (file_exists) {
            std::getline(check, first_line);
        }
        EXPECT_EQ(first_line, "TITL PDB->ins");
    }

    // Case 2: mol_is_from_shelx_ins = true
    {
        std::string outfile = "/tmp/test_shelx_ins_2.ins";
        coot::ShelxIns shelx;
        std::pair<int, std::string> result = shelx.write_ins_file_gemmi(st, outfile, true);

        EXPECT_EQ(result.first, 1);
        EXPECT_EQ(result.second, "INFO:: SHELXL file /tmp/test_shelx_ins_2.ins written.");

        std::ifstream check(outfile);
        bool file_exists = check.good();
        EXPECT_TRUE(file_exists);

        std::string first_line;
        if (file_exists) {
            std::getline(check, first_line);
        }
        EXPECT_EQ(first_line, "TITL PDB->ins");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}