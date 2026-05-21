#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, copy_cell_and_symm_headers) {
    // case: copy_cell_and_symm_headers_gemmi copies crystallographic data from one Structure to another
    
    gemmi::Structure src = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Structure dest = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Get initial crystallographic data from dest
    double cell_a = dest.cell.a;
    double cell_b = dest.cell.b;
    double cell_c = dest.cell.c;
    double cell_alpha = dest.cell.alpha;
    double cell_beta = dest.cell.beta;
    double cell_gamma = dest.cell.gamma;
    const std::string& sg_before = dest.spacegroup_hm;
    
    EXPECT_NEAR(cell_a, 53.512, 1e-4);
    EXPECT_NEAR(cell_b, 71.311, 1e-4);
    EXPECT_NEAR(cell_c, 72.101, 1e-4);
    EXPECT_NEAR(cell_alpha, 90.0, 1e-4);
    EXPECT_NEAR(cell_beta, 90.0, 1e-4);
    EXPECT_NEAR(cell_gamma, 90.0, 1e-4);
    EXPECT_STREQ(sg_before.c_str(), "P 21 21 21");
    
    // Copy crystallographic data from src to dest
    bool result = coot::util::copy_cell_and_symm_headers_gemmi(src, dest);
    
    EXPECT_TRUE(result);
    
    // Get crystallographic data from dest after copy
    cell_a = dest.cell.a;
    cell_b = dest.cell.b;
    cell_c = dest.cell.c;
    cell_alpha = dest.cell.alpha;
    cell_beta = dest.cell.beta;
    cell_gamma = dest.cell.gamma;
    const std::string& sg_after = dest.spacegroup_hm;
    
    EXPECT_NEAR(cell_a, 53.512, 1e-4);
    EXPECT_NEAR(cell_b, 71.311, 1e-4);
    EXPECT_NEAR(cell_c, 72.101, 1e-4);
    EXPECT_NEAR(cell_alpha, 90.0, 1e-4);
    EXPECT_NEAR(cell_beta, 90.0, 1e-4);
    EXPECT_NEAR(cell_gamma, 90.0, 1e-4);
    EXPECT_STREQ(sg_after.c_str(), "P 21 21 21");
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}