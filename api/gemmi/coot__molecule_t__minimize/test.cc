#include <gtest/gtest.h>
#include <string>
#include "function.hh"

TEST(OracleTest, MoleculeMinimizeStatus) {
   coot::molecule_t mol;
   int status = mol.minimize_gemmi("//A/10", 
                                  5, 
                                  false, 0.5, 
                                  false, 0.5, 
                                  true,
                                  nullptr);
   EXPECT_EQ(status, -2);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}