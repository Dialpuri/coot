#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_group) {
    coot::protein_geometry geom;
    geom.init_standard();

    // Case 1: standard amino acid "ALA"
    {
        std::string res_name = "ALA";
        std::string group = coot::get_group_gemmi(geom, res_name);
        EXPECT_EQ(group, "peptide");
    }

    // Case 2: nucleic acid residue "A" (adenine)
    {
        std::string res_name = "A";
        std::string group = coot::get_group_gemmi(geom, res_name);
        EXPECT_EQ(group, "RNA");
    }

    // Case 3: invalid residue name — should throw
    {
        std::string res_name = "ZZZ";
        EXPECT_THROW(coot::get_group_gemmi(geom, res_name), std::runtime_error);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
