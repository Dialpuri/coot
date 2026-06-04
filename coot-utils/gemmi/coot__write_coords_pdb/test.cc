#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include <fstream>

TEST(OracleTest, write_coords_pdb) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_GT(st.models.size(), 0u);  // ensure it loaded

    // Case 1: valid molecule, write to a temp file
    {
        std::string file_name = "/tmp/oracle_write_coords_pdb_case1.pdb";

        int result = coot::write_coords_pdb_gemmi(st, file_name);

        std::ifstream ifs(file_name);
        bool file_exists = ifs.good();
        ifs.close();

        EXPECT_EQ(result, 0);
        EXPECT_TRUE(file_exists);
    }

    // Case 2: write to a different temp file to confirm it works again
    {
        std::string file_name = "/tmp/oracle_write_coords_pdb_case2.pdb";

        int result = coot::write_coords_pdb_gemmi(st, file_name);

        std::ifstream ifs(file_name);
        bool file_exists = ifs.good();
        ifs.close();

        EXPECT_EQ(result, 0);
        EXPECT_TRUE(file_exists);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}