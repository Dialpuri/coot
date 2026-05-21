#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, copy_headers) {
    std::string pdb_path = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";
    
    gemmi::Structure st1 = gemmi::read_pdb_file(pdb_path);
    gemmi::Structure st2 = gemmi::read_pdb_file(pdb_path);
    
    int nRemarks1 = st1.raw_remarks.size();
    int nRemarks2_before = st2.raw_remarks.size();
    
    EXPECT_EQ(nRemarks1, 94);
    EXPECT_EQ(nRemarks2_before, 94);
    
    {
        bool include_cryst = false;
        bool result = coot::util::copy_headers_gemmi(&st1, &st2, include_cryst);
        EXPECT_TRUE(result);
        
        int nRemarks2_after = st2.raw_remarks.size();
        EXPECT_EQ(nRemarks2_after, 94);
    }
    
    st2 = gemmi::read_pdb_file(pdb_path);
    
    {
        bool include_cryst = true;
        bool result = coot::util::copy_headers_gemmi(&st1, &st2, include_cryst);
        EXPECT_TRUE(result);
        
        int nRemarks2_after = st2.raw_remarks.size();
        EXPECT_EQ(nRemarks2_after, 94);
    }
    
    {
        bool include_cryst = false;
        bool result = coot::util::copy_headers_gemmi(nullptr, &st2, include_cryst);
        EXPECT_FALSE(result);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}