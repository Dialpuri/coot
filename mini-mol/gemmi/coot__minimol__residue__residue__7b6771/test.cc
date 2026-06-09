#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, residue) {
    // Case 1: construct with seqnum and resname
    {
        int seqnum_input = 42;
        std::string resname_input = "ALA";
        coot::minimol::residue r = coot::minimol::residue_gemmi(seqnum_input, resname_input);
        EXPECT_EQ(r.seqnum, 42);
        EXPECT_EQ(r.name, "ALA");
        EXPECT_EQ(r.ins_code, "");
        EXPECT_EQ(r.n_atoms(), 0u);
    }

    // Case 2: construct with a different residue (GLY, seqnum 1)
    {
        int seqnum_input = 1;
        std::string resname_input = "GLY";
        coot::minimol::residue r = coot::minimol::residue_gemmi(seqnum_input, resname_input);
        EXPECT_EQ(r.seqnum, 1);
        EXPECT_EQ(r.name, "GLY");
        EXPECT_EQ(r.ins_code, "");
        EXPECT_EQ(r.n_atoms(), 0u);
    }

    // Case 3: construct with seqnum 0 and empty resname (edge case)
    {
        int seqnum_input = 0;
        std::string resname_input = "";
        coot::minimol::residue r = coot::minimol::residue_gemmi(seqnum_input, resname_input);
        EXPECT_EQ(r.seqnum, 0);
        EXPECT_EQ(r.name, "");
        EXPECT_EQ(r.ins_code, "");
        EXPECT_EQ(r.n_atoms(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
