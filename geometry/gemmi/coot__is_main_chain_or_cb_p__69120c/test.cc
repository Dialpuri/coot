#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, is_main_chain_or_cb_p) {
    // Case 1: valid main chain atom (CA)
    {
        std::string atom_name = " CA ";
        bool result = coot::is_main_chain_or_cb_p_gemmi(atom_name);
        EXPECT_TRUE(result);
    }

    // Case 2: CB atom (should also be true)
    {
        std::string atom_name = " CB ";
        bool result = coot::is_main_chain_or_cb_p_gemmi(atom_name);
        EXPECT_TRUE(result);
    }

    // Case 3: side chain atom not in the list (should be false)
    {
        std::string atom_name = " CG ";
        bool result = coot::is_main_chain_or_cb_p_gemmi(atom_name);
        EXPECT_FALSE(result);
    }

    // Case 4: nitrogen (main chain)
    {
        std::string atom_name = " N  ";
        bool result = coot::is_main_chain_or_cb_p_gemmi(atom_name);
        EXPECT_TRUE(result);
    }

    // Case 5: invalid/unrecognized atom
    {
        std::string atom_name = " OH ";
        bool result = coot::is_main_chain_or_cb_p_gemmi(atom_name);
        EXPECT_FALSE(result);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}